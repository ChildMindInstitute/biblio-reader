import pandas as pd
import unidecode, os, numpy
import xml.etree.ElementTree as etree
from Bio import Entrez
import manager
Entrez.email = 'drcc@vt.edu'
pd.options.mode.chained_assignment = None


def filterstr(str, filter, decode=True):
    for char in filter:
        str = str.replace(char, '')
    if decode:
        return unidecode.unidecode(str)
    else:
        return str


def get_ids(data):
    ids = []
    for row in data.iterrows():
        row = row[1]
        year = str(row['Year']).replace('.0', '')
        term = filterstr(str(row['Title']), ")'>',\```[](<}{")
        author = filterstr(str(row['Authors']).lower().split(' & ')[0], ")'',```'[\](}{')")
        request = Entrez.esearch(db='pubmed', term=term, retmax=1, field='title', year=year, author=author)
        idlist = Entrez.read(request)['IdList']
        if len(idlist) == 1:
            ids.append(int(idlist[0]))
        else:
            ids.append(0)
    data['PMCID'] = ids
    return ids


def write_bib(data, directory):
    if 'PMCID' not in data:
        print("NO PMCIDS")
        return
    id_dict = {i: int(pmcid) for i, pmcid in dict(zip(data['i'], data['PMCID'])).items() if int(pmcid) != 0}
    for i, pmcid in id_dict.items():
        bib = Entrez.efetch(db='pubmed', id=pmcid, retmode="xml", rettype="full")
        with open('/'.join([directory, str(i) + '.xml']), 'w') as f:
            f.write(bib.read())


def parse_bib(directory):
    bibs = {}
    for bib in os.listdir(directory):
        root = etree.parse(open('/'.join([directory, bib]))).getroot()
        authors = ''
        for auth in root.findall('.//Author'):
            fore = auth.find('ForeName')
            last = auth.find('LastName')
            if fore is not None and last is not None:
                authors += ' '.join([fore.text, last.text, '&'])
        affiliations = ' && '.join(set([aff.text for aff in root.findall('.//Affiliation')]))
        qualifiers = ' & '.join(set([qual.text for qual in root.findall('.//MeshHeading/QualifierName')]).union(set([
            key.text for key in root.findall('.//KeywordList/Keyword')])))
        bibs[int(bib.replace('.xml', ''))] = (authors, affiliations, qualifiers)
    bibs.update({i: (None, numpy.NaN, numpy.NaN) for i in range(0, 1560) if i not in bibs.keys()})
    return {key: bib for key, bib in sorted(bibs.items())}




if __name__ == '__main__':
    get_ids(manager.DATA)
    manager.update()
    write_bib(manager.DATA, manager)