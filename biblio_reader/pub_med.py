import pandas as pd, unidecode, os, xml.etree.ElementTree as etree, sys
# from biopython <http://biopython.org/DIST/docs/api/Bio.Entrez-module.html>
from Bio import Entrez
br_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
if br_path not in sys.path:
    sys.path.append(br_path)
Entrez.email = 'drcc@vt.edu'
pd.options.mode.chained_assignment = None
import manager as mg
import urllib.request as urllib


def filterstr(str, filter, decode=True):
    """
    Returns a string without characters contained in the filter
    :param str: The string to be filtered
    :param filter: The characters that will filter the string
    :param decode: If true, returns a UTF-8 encoded string
    :return: Filtered string
    """
    for char in filter:
        str = str.replace(char, '')
    if decode:
        return unidecode.unidecode(str)
    else:
        return str
    

def get_abstract(data):
    """
    Function to get a text abstract from any article with a PubMed ID
    
    Parameter
    ---------
    pmid : int
    
    Returns
    -------
    abstract : string
    """

    abstracts = []
    pmcids = list(data['PMCID'])

    for int in range(len(pmcids)):
        pmid = pmcids[int]

        try:
            abstract = "Not in PubMed" if pmid == 0 else urllib.urlopen(
                "http://togows.dbcls.jp/entry/ncbi-pubmed/{0}/abstract".format(
                pmid)
            ).read().decode("UTF-8")
        except:
            abstract = "No abstract"

        abstracts.append(abstract)

    data['Abstracts'] = abstracts
    return(abstracts)

def get_journals(data):

    journals = []
    pmcids = list(data['PMCID'])

    for int in range(len(pmcids)):
        pmid = pmcids[int]

        try:
            journal = "Not in PubMed" if pmid == 0 else urllib.urlopen(
                "http://togows.dbcls.jp/entry/ncbi-pubmed/{0}/journal".format(
                    pmid)
            ).read().decode("UTF-8")
        except:
            journal = "No Journal Found"

        journals.append(journal)

    data['Journal'] = journals

    return journals

def get_ids(data):
    """
    Takes all pubs of the dataframe and searches Pubmed for PMCID's related to
    the title, authors, and year in each row
    :param data: The dataframe
    :return: a list of all IDS, an updated dataframe
    """
    ids = []
    for row in data.iterrows():
        row = row[1]
        year = row['Year'] or ''
        term = filterstr(str(row['Title']), ")'>',\```[](<}{")
        author = filterstr(str(row['Authors']).lower().split(' & ')[0],
                 ")'',```'[\](}{')")
        request = Entrez.esearch(db='pubmed', term=term, retmax=1, field=
                  'title', year=year, author=author)
        idlist = Entrez.read(request)['IdList']
        if len(idlist) == 1:
            ids.append(int(idlist[0]))
        else:
            ids.append(0)
    data['PMCID'] = ids
    return ids


def write_bib(data, directory):
    """
    Takes all PMCIDS and gets the associated bibliographies for each
    :param data: The dataframe containing PMCIDS
    :param directory: The directory to write bibs in
    """
    if 'PMCID' not in data:
        print("NO PMCIDS")
        return
    id_dict = {i: int(pmcid) for i, pmcid in dict(zip(data['i'], data['PMCID'])
              ).items() if int(pmcid) != 0}
    for i, pmcid in id_dict.items():
        bib = Entrez.efetch(db='pubmed', id=pmcid, retmode="xml", rettype=
              "full")
        with open('/'.join([directory, str(i) + '.xml']), 'w') as f:
            f.write(bib.read())


def parse_bib(directory, outfile):
    """
    Parses each Pubmed obtained bibliography to find qualifiers, affiliations,
    and full author lists
    :param directory: The directory of Pubmed bibs
    :param outfile: outfile dir to write result in
    :return: A CSV file with reference #s, qualifiers, affiliations, and
             authors for each article with a PMCID
    """
    parsed = []
    for bib in os.listdir(directory):
        root = etree.parse(open('/'.join([directory, bib]))).getroot()
        authors = []
        for auth in root.findall('.//Author'):
            fore = auth.find('ForeName')
            last = auth.find('LastName')
            if fore is not None and last is not None:
                authors.append(' '.join([fore.text, last.text]))
        authors = ';;'.join(authors)
        affiliations = ';;'.join(set([aff.text for aff in root.findall(
                       './/Affiliation')]))
        qualifiers = ';;'.join(set([qual.text for qual in root.findall(
                     './/MeshHeading/QualifierName')]).union(set([
            key.text for key in root.findall('.//KeywordList/Keyword')])))
        parsed.append((int(bib.replace('.xml', '')), authors, affiliations,
                      qualifiers))
    parsed_data = pd.DataFrame(parsed, columns=['i', 'authors', 'affiliations',
                  'qualifiers'])
    parsed_data.sort_values('i', inplace=True)
    parsed_data.to_csv(path_or_buf=outfile, index=False)


if __name__ == '__main__':

    data = mg.get_data()
    BIB_DIR = mg.dir(os.path.join(mg.WORKING_PATH, 'bibs'))
    PARSED_BIBS = os.path.join(mg.WORKING_PATH, 'parsed_bibs.csv')

    if 'PMCID' not in data:
        get_ids(data)
        mg.update_data()
    if not os.path.exists(BIB_DIR):
        write_bib(data, mg.dir(BIB_DIR))
    parse_bib(BIB_DIR, PARSED_BIBS)
    # get_abstract(data)
    get_journals(data)
    mg.update_data()
