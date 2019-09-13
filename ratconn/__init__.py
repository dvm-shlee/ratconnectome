import os
import re
import uncurl
import requests
from pandas import DataFrame
from bs4 import BeautifulSoup

__version__ = '0.0.2'

#%%
curl_cmd = "curl 'http://neuroviisas.med.uni-rostock.de/connectome/index.php' \
-XPOST \
-H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8' \
-H 'Content-Type: application/x-www-form-urlencoded' \
-H 'Origin: http://neuroviisas.med.uni-rostock.de' \
-H 'Content-Length: 94' \
-H 'Accept-Language: en-us' \
-H 'Upgrade-Insecure-Requests: 1' \
-H 'Host: neuroviisas.med.uni-rostock.de' \
-H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_6) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/12.1.2 Safari/605.1.15' \
-H 'Referer: http://neuroviisas.med.uni-rostock.de/connectome/index.php' \
-H 'Accept-Encoding: gzip, deflate' \
-H 'Connection: keep-alive' \
--data 'searchOpt=starts&searchSide=searchAll&searchNameText=&searchName=Search+name&searchAbbrevText='"


def parse_request_cmd(curl_cmd):
    request_cmd = uncurl.parse(curl_cmd)
    request_text = "".join(l.strip() for l in request_cmd.split('\n'))
    p_url = re.compile(r'requests.post\(((?!,).)*', )
    p_headers = re.compile(r'headers={((?!{).)*}')
    output = dict()

    u_start, u_end = p_url.search(request_text).span()
    h_start, h_end = p_headers.search(request_text).span()

    output["url"] = request_text[u_start + 15:u_end - 1]
    exec(request_text[h_start:h_end], output)
    return output['headers'], output['url']


def parse_tables_content(html):
    soup = BeautifulSoup(html, 'html.parser')
    tables = soup.find_all('table', attrs={'class':'tableStyle'})
    if len(tables) > 1:
        return tables
    else:
        return tables[0]


def parse_table(table):
    has_url = True if len(table.find_all('a', attrs={'href': re.compile("^showRegion")})) != 0 else False
    rows = table.find_all('tr')
    output, columns = None, None
    for i, r in enumerate(rows):
        if i is 0:
            columns = [th.text.strip('.') for th in r.find_all('th')]
            if has_url:
                columns += ['Url']
            output = DataFrame(columns = columns)
        else:
            values = list()
            for td in r.find_all('td'):
                if len(td.text) > 0:
                    values.append(td.text)
                else:
                    values.append(td.get('bgcolor'))
            if has_url:
                values += [r.find('a', attrs={'href': re.compile("^showRegion")}).get('href')]
            output = output.append(dict(zip(columns, values)), ignore_index=True)
    return output


def get_conn_matrix(*rois):
    itf = Interface(abbr=True, verbose=False)
    columns = ['{}-left'.format(roi) for roi in rois]
    columns += ['{}-right'.format(roi) for roi in rois]
    output_df = DataFrame(index=columns, columns=columns)
    input_df = DataFrame(index=columns, columns=columns)
    sides = ['left', 'right']
    for roi in rois:
        for i, s in enumerate(sides):
            itf.search(roi)
            row_item = '{}-{}'.format(roi, s)
            info = itf.get_info(i)
            ### output_matrix
            o_df = info.output_to
            o_df = o_df.loc[o_df['Abbr'].isin(list(rois))].reset_index(drop=True)
            i_df = info.input_from
            i_df = i_df.loc[i_df['Abbr'].isin(list(rois))].reset_index(drop=True)
            for j, row in o_df.iterrows():
                if row.Side == 'ipsi':
                    side = sides[i]
                else:
                    side = sides[i-1]
                col_item = '{}-{}'.format(row.Abbr, side)
                output_df.loc[row_item, col_item] = float(row.Weight)
            for j, row in i_df.iterrows():
                if row.Side == 'ipsi':
                    side = sides[i]
                else:
                    side = sides[i-1]
                col_item = '{}-{}'.format(row.Abbr, side)
                input_df.loc[row_item, col_item] = float(row.Weight)
    result = output_df.where(output_df > input_df, input_df).fillna(output_df)
    return result.fillna(0)


class Info():
    def __init__(self, url, *args, **kwargs):
        self.verbose = kwargs['verbose'] \
            if 'verbose' in kwargs.keys() and isinstance(kwargs['verbose'], bool) else True
        self._base_url = url
        self.synonyms = args[0]
        self.subregions = args[1]
        self._output_to = args[2]
        self._input_from = args[3]
        self.weight = args[4]
        self.source = 'output'

    @property
    def output_to(self):
        output_df = self._apply_filter(self._output_to)
        return output_df.loc[~output_df['Weight'].isin(['0'])].reset_index(drop=True)

    @property
    def input_from(self):
        input_df = self._apply_filter(self._input_from)
        return input_df.loc[~input_df['Weight'].isin(['0'])].reset_index(drop=True)

    def _apply_filter(self, df):
        df = df.sort_values(by=['Weight', 'Name'], ascending=False)
        return df

    def get_info(self, idx):
        if self.source == 'output':
            table = self.output_to
        else:
            table = self.input_from
        url = os.path.join(self._base_url, table.loc[idx, 'Url'])
        html = requests.get(url).content
        soup = BeautifulSoup(html, 'html.parser')
        tables = soup.find_all('table', attrs={'class': 'tableStyle'})
        output = []
        for i, t in enumerate(tables):
            if i != 0:
                output.append(parse_table(t))

        return Info(self._base_url, *output)

    def toggle_source(self):
        if self.source == 'output':
            self.source = 'input'
        else:
            self.source = 'output'

class Interface():
    def __init__(self, roi=None, **kwargs):
        sides = ['searchAll', 'left', 'right']
        options = ['starts', 'contains']

        self.headers, self.url = parse_request_cmd(curl_cmd)
        self._base_url = os.path.dirname(self.url)
        self.side = kwargs['side'].lower() \
            if ('side' in kwargs.keys() and kwargs['side'] in sides) else 'searchAll'
        self.option = kwargs['option'].lower() \
            if ('option' in kwargs.keys()) and kwargs['option'] in options else 'starts'
        self.abbr = kwargs['abbr'] \
            if 'abbr' in kwargs.keys() and isinstance(kwargs['abbr'], bool) else  False
        self.verbose = kwargs['verbose'] \
            if 'verbose' in kwargs.keys() and isinstance(kwargs['verbose'], bool) else True
        if roi is not None:
            self.search(roi)
        else:
            self.roi = None

    def search(self, roi):
        self.roi = roi
        html = requests.post(url=self.url, data=self._get_data_dict()).content
        table_content = parse_tables_content(html)
        self.results = parse_table(table_content)

    def get_info(self, idx):
        url = os.path.join(self._base_url, self.results.loc[idx, 'Url'])
        html = requests.get(url).content
        soup = BeautifulSoup(html, 'html.parser')
        tables = soup.find_all('table', attrs={'class': 'tableStyle'})
        output = []
        for i, t in enumerate(tables):
            if i != 0:
                output.append(parse_table(t))
        if self.verbose:
            print("Connectivity information of {}".format(self.results.loc[idx, 'Name']))
        return Info(self._base_url, *output)

    def toggle_side(self):
        if self.side == 'searchAll':
            self.side = 'left'
        elif self.side == 'left':
            self.side = 'right'
        elif self.side == 'right':
            self.side = 'searchAll'
        else:
            pass
        if self.verbose:
            print('Toggled to search "{}" side.'.format(self.side))

    def toggle_option(self):
        if self.option == 'starts':
            self.option = 'contains'
            if self.verbose:
                print('Search the ROI contains "{}".'.format(self.roi))
        else:
            self.option = 'starts'
            if self.verbose:
                print('Search the ROI start with "{}".'.format(self.roi))

    def toggle_abbr(self):
        if self.abbr:
            self.abbr = False
            if self.verbose:
                print('Search with full name.')
        else:
            self.abbr = True
            if self.verbose:
                print('Search with abbreviation.')

    def _get_data_dict(self):
        output = dict(searchOpt=self.option,
                      searchSide=self.side)
        if self.abbr is True:
            output['searchNameText'] = None
            output['searchAbbrevText'] = self.roi
            output['searchAbbrev'] = 'Search+abbreviation'
        else:
            output['searchNameText'] = self.roi
            output['searchAbbrevText'] = None
            output['searchName'] = 'Search+name'
        return output

    def __repr__(self):
        line = ['Input ROI = {}'.format(self.roi)]
        line.append('Side = {}'.format(self.side))
        line.append('Search option = {}'.format(self.option))
        line.append('Abbr. = {}'.format(self.abbr))
        return '\n'.join(line)