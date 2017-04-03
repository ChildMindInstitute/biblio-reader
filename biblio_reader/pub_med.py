import pandas as pd
import unidecode
import urllib.request as urllib
import urllib.error as urlerror
import xml.etree.ElementTree as etree
import re
import webbrowser
from bs4 import BeautifulSoup
pd.options.mode.chained_assignment = None

with open('inputs/FCP_DATA_wIDS.csv', 'r') as f:
    data = pd.read_csv(f)


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
    res = []
    try:
        xml = urllib.urlopen(url).read()
    except:
        return res
    for idlist in etree.fromstring(xml).findall('IdList'):
        for id in idlist:
            res.append(int(id.text))
    return res


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
    keywords = map(lambda x: re.compile(x), keywords)
    for link in fulldata['Full']:
        try:
            html = urllib.urlopen(link).read()
        except urlerror:
            continue
        soup = BeautifulSoup(html, 'html.parser')
        excerpts.append(soup.find_all(string=keywords))
    fulldata['Excerpts'] = excerpts





data['PMCIDS'] = data['PMCIDS'].apply(change_id)
data['PMCIDS'] = data['PMCIDS'].fillna(0.0).astype(int)

fulldata = data[data['PMCIDS'] != 0]

fulldata.reset_index(inplace=True)
del fulldata['Unnamed: 0']

fulldata['Full'] = fulldata['PMCIDS'].apply(lambda x: 'https://eutils.ncbi.nlm.nih' +
                                                      '.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=' +
                                                      str(x) + '&cmd=prlinks&retmode=ref')
with open('outputs/FCP_DATA_wIDS.csv', 'w') as f:
    f.write(fulldata.to_csv(index=False))


fcp = """"fcon_1000.projects.nitrc.org" OR "Rockland Sample" OR "1000 Functional Connectomes Project" OR "International Neuroimaging Data-Sharing Initiative" OR "Autism Brain Imaging Data Exchange" OR ADHD-200 OR "Consortium for Reproducibility and Reliability”""".replace('"', '')

get_key_paragraphs(fcp.split(' OR '))
