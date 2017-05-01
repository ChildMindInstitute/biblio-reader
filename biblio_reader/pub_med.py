import pandas as pd
import unidecode
import urllib.request as urllib
import xml.etree.ElementTree as etree
import re
from bs4 import BeautifulSoup
pd.options.mode.chained_assignment = None

with open('../inputs/FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)


fcp = ['fcon_1000.projects.nitrc.org', 'Rockland Sample', '1000 Functional Connectomes',
       'International Neuroimaging Data-Sharing Initiative', 'Autism Brain Imaging Data Exchange', 'ADHD-200',
       'Consortium for Reproducibility and Reliability']


def safeint(x):
    try:
        return int(x)
    except:
        return x


def filterstr(str, filter, decode=True):
    for char in str:
        if char in filter:
            str = str.replace(char, '')
    if decode:
        return unidecode.unidecode(str)
    else:
        return str


def joindata(row):
    row.fillna('')
    return '&'.join(['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed', 'term='
                     + filterstr(str(row['Title']).strip().lower(), ")'>',\```[](<}{"), 'retmax=1', 'field=title',
                     'year=' + str(safeint(row['Year'])), 'author=' +
                     filterstr(str(row['Authors']).lower().split(' & ')[0],
                               ")'',```'[\](}{')")]).replace(' ', '+').replace('"', '')


def get_id(url):
    try:
        xml = urllib.urlopen(url).read()
    except:
        return None
    for idlist in etree.fromstring(xml).findall('IdList'):
        for id in idlist:
            return int(id.text)


def change_id(id):
    if isinstance(id, list):
        if len(id) == 1:
            return id[0]
    return id


def get_ids(file=None, full_text=True):
    ids = []
    for row in data.iterrows():
        ids.append(getid(joindata(row[1])))
    data['PMCIDS'] = ids
    data['PMCIDS'] = data['PMCIDS'].apply(change_id).fillna(0.0).astype(int)
    fulldata = data[data['PMCIDS'] != 0]
    fulldata['index'] = list(range(0, len(fulldata['index'])))
    fulldata.reset_index(inplace=True)
    if full_text:
        fulldata['Full'] = fulldata['PMCIDS'].apply(lambda x: 'https://eutils.ncbi.nlm.nih' +
                                                              '.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=' +
                                                              str(x) + '&cmd=prlinks&retmode=ref')
    if file:
        fulldata.to_csv(path_or_buf=file)


def get_key_paragraphs(keywords):
    if fulldata['Full'] is None:
        return
    excerpts = []
    keywords = re.compile("(" + '|'.join(keywords) + ")")
    for link in fulldata['Full']:
        try:
            html = urllib.urlopen(link).read()
        except:
            excerpts.append(['N/A'])
            continue
        soup = BeautifulSoup(html, 'html.parser')
        excerpts.append(soup.find_all(string=keywords))
    fulldata['Excerpts'] = excerpts

fulldata['index'] = list(range(0, len(fulldata['index'])))


with open('outputs/FCP_DATA_wIDS.csv', 'w') as f:
    f.write(fulldata.to_csv(index=False))

