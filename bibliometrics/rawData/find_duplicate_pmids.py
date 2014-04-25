from Bio import Entrez
import re
from pybtex.database.input import bibtex
from sets import Set

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

pubmed_ids = Set()
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
    else:
        if bib_data.entries[k].fields['pmid'] not in pubmed_ids:
             pubmed_ids.add(bib_data.entries[k].fields['pmid'])
        else:
            print "%s is a duplicate"%(bib_data.entries[k].fields['pmid'])

print "%d unique pubmed ids"%(len(pubmed_ids))
