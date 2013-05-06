from Bio import Entrez
import sets
from pybtex.database.input import bibtex

def is_int(a):
    try:
        int (a)
        return True
    except:
        return False

## setup access to pubmed
Entrez.email='drcc@vt.edu'
database='pubmed'

## setup bibtex parser
parser = bibtex.Parser()

bib_data=parser.parse_file('cmi_library.bib')

id_num=0
lib_pmids=set()

for k in bib_data.entries.keys():
    if 'pmid' in bib_data.entries[k].fields.keys():
        lib_pmids.add(bib_data.entries[k].fields['pmid'])
    else:
        print "PubMed ID missing for %s: %s"%(k, \
            bib_data.entries[k].fields['title'])

# now we perform the new query
query = """
    ("rest" or "resting state" or "spontaneous" or "intrinsic" or
     "low frequency fluctuation" or "low frequency fluctuations" or 
     "connectivity" or "connectome" or "default network" or 
     "default mode network" or "default network") and ("functional magnetic resonance 
     imaging" or "functional mri" or "fmri") and ("1901/01/01"[Date - 
     Publication] : "2012/12/31"[Date - Publication])
"""

print query


result = Entrez.esearch(db=database,term=(query),retmax=100000)
try:
    record=Entrez.read(result)
except:
    print "Error"

print "\n\nSearch returned %d records"%(int(record['Count']))

new_pmids=set(record['IdList'])

new_diff_pmids=new_pmids.difference(lib_pmids)
old_diff_pmids=lib_pmids.difference(new_pmids)

print "Received %d new pmids"%(len(new_pmids))
with open('new_pmids.txt', 'w') as f:
    f.write("\n".join(new_pmids))

print "%d of which are not in the old library"%(len(new_diff_pmids))
with open('new_pmids_diff.txt', 'w') as f:
    f.write("\n".join(new_diff_pmids))

print "%d pmids in the old library are not in the new list"%(len(old_diff_pmids))
with open('old_pmids_diff.txt', 'w') as f:
    f.write("\n".join(old_diff_pmids))

