from Bio import Entrez
import re
from pybtex.database.input import bibtex
#from pybtex.database.output.bibtex import Writer

#def is_int(a):
    #try:
        #int (a)
        #return True
    #except:
        #return False

## setup access to pubmed
#Entrez.email='drcc@vt.edu'
#database='pubmed'

## setup bibtex parser
parser = bibtex.Parser()
## setup bibtex write
#w = Writer(); 

bib_data=parser.parse_file('cmi_library.bib')

id_num=0
for k in bib_data.entries.keys():
    if not 'pmid' in bib_data.entries[k].fields.keys():
        title="%s"% \
            (re.sub('[\{\}\\\\]','',bib_data.entries[k].fields['title']))
        print "Key %s is missing: %s"%(k,title)
