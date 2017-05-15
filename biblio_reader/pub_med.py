import pandas as pd
import unidecode
import xml.etree.ElementTree as etree
import os
from Bio import Entrez
Entrez.email = 'drcc@vt.edu'
pd.options.mode.chained_assignment = None

with open('../inputs/FCP_DATA.csv', 'r') as f, open('../outputs/PMCIDS.txt', 'r') as g:
    data = pd.read_csv(f)
    pmcids = pd.read_csv(g)


def safeint(x):
    try:
        return int(x)
    except:
        return x


def filterstr(str, filter, decode=True):
    for char in filter:
        str = str.replace(char, '')
    if decode:
        return unidecode.unidecode(str)
    else:
        return str


def get_ids(file=None):
    ids = []
    for row in data.iterrows():
        row = row[1]
        year = str(safeint(row['Year']))
        term = filterstr(str(row['Title']), ")'>',\```[](<}{")
        author = filterstr(str(row['Authors']).lower().split(' & ')[0], ")'',```'[\](}{')")
        request = Entrez.esearch(db='pubmed', term=term, retmax=1, field='title', year=year, author=author)
        idlist = Entrez.read(request)['IdList']
        if len(idlist) == 1:
            ids.append(int(idlist[0]))
        else:
            ids.append(0)
    data['PMCIDS'] = ids
    if file:
        data.to_csv('../inputs/FCP_DATA.csv')
    return ids


def write_bib(directory):
    id_dict = {i: int(pmcid) for i, pmcid in dict(zip(data['i'], data['PMCIDS'])).items() if int(pmcid) != 0}
    for i, pmcid in id_dict.items():
        bib = Entrez.efetch(db='pubmed', id=pmcid, retmode="xml", rettype="full")
        with open('/'.join([directory, str(i) + '.xml']), 'w') as f:
            f.write(bib.read())


for bib in os.listdir('../outputs/bibs'):
    root = etree.parse(open('../outputs/bibs/' + bib)).getroot()
    title = root.findall('.//ArticleTitle')[0].text.lower()
    title2 = data.iloc[int(bib.replace('.xml', ''))]['Title'].lower() + '.'
    if title2 != title:
        print(bib.replace('.xml', ''), title, title2)
