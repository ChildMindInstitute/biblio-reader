from Bio import Entrez
import sets
from pybtex.database.input import bibtex

## setup access to pubmed
Entrez.email='drcc@vt.edu'
database='pubmed'

def getNames(authorList):
    names = []
    for author in authorList:
        if 'ForeName' not in author:
            if 'LastName' in author:
                names.append(author['LastName'])
            else:
                assert 'CollectiveName' in author
                names.append(author['CollectiveName'])

        elif 'LastName' not in author:
            assert 'ForeName' in author
            names.append(author['ForeName'])
        else:
            names.append(author['LastName'] + ", " + author['ForeName'])
    ret = " and ".join(names)
    return ret

def writeBib(inFile, outFile):
    out = open(outFile, 'w')
    for line in open(inFile, 'r').readlines():
        pmid = line.rstrip()
        k_result = Entrez.efetch(db=database, id=pmid, retmode="xml", \
                                 rettype="full", tool="bioPython, Bio.Entrez")
        record = Entrez.read(k_result)
        if len(record) != 1 or 'MedlineCitation' not in record[0] or 'Article' not in record[0]['MedlineCitation']:
            print "MISSING PUBMED DATA FOR", pmid
            continue

        e = {}
        e['author'] = getNames(record[0]['MedlineCitation']['Article']['AuthorList'])
        e['title'] = record[0]['MedlineCitation']['Article']['ArticleTitle']
        if 'Abstract' in record[0]['MedlineCitation']['Article']:
            e['abstract'] = record[0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
        if 'ELocationID' in record[0]['MedlineCitation']['Article'] and len(record[0]['MedlineCitation']['Article']['ELocationID']) > 0:
            e['doi'] = record[0]['MedlineCitation']['Article']['ELocationID'][0]
        if 'ISSN' in record[0]['MedlineCitation']['Article']['Journal']:
            e['isbn'] = record[0]['MedlineCitation']['Article']['Journal']['ISSN']
        e['language'] = ""
        e['url'] = ""
        if 'Issue' in record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']:
            e['number'] = record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['Issue']
        e['pages'] = record[0]['MedlineCitation']['Article']['Pagination']['MedlinePgn']
        if 'Volume' in record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']:
            e['volume'] = record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['Volume']
        e['file'] = ""
        if 'Year' in record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            e['year'] = record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
        elif 'MedlineDate' in record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
            e['year'] = str(int(record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'].split(" ")[0][:4]))
        else:
            print "NO YEAR", record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
        if 'MeshHeadingList' in record[0]['MedlineCitation']:
            e['mesh_heading_list'] = " and ".join([x['DescriptorName'] for x in record[0]['MedlineCitation']['MeshHeadingList']])

        e['journal'] = record[0]['MedlineCitation']['Article']['Journal']['Title']
        e['journal_iso'] = record[0]['MedlineCitation']['Article']['Journal']['ISOAbbreviation']
        if 'Affiliation' in record[0]['MedlineCitation']['Article']:
            e['affiliation'] = record[0]['MedlineCitation']['Article']['Affiliation']
        
        firstLastName = record[0]['MedlineCitation']['Article']['AuthorList'][0]['LastName']
        bibId = firstLastName + e['year']
        out.write("@article{\n")
        out.write(bibId.encode('utf8') + ",\n")
        entries = []
        for key in e.keys():
            if len(e[key]) > 0:
                if isinstance(e[key], str) or isinstance(e[key], unicode):
                    x = e[key].encode('utf8')
                elif isinstance(e[key], Entrez.Parser.StringElement):
                    x = str(e[key])
                elif isinstance(e[key], Entrez.Parser.StructureElement):
                   x = str(e[key])
                elif isinstance(e[key], list):
                    x = str(e[key])   
                else:
                    print key
                    print "Type is " + str(type(e[key]))
                    print str(e[key])
                    assert False
                entries.append(key + ' = "' + x + '"')
        out.write(",\n".join(entries))
        out.write("\n}\n")
    pass

if __name__ == "__main__":
    for f in ["new_pmids_diff.txt", "new_pmids.txt"]:#,"old_pmids_diff.txt"]:
        writeBib(f, f.split(".")[0] + ".bib")
    