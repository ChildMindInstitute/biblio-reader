from Bio import Entrez
import re
from pybtex.database.input import bibtex
from pybtex.database.output.bibtex import Writer

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
## setup bibtex write
w = Writer(); 

bib_data=parser.parse_file('cmi_library.bib')

id_num=0
for k in bib_data.entries.keys():
    if not 'pmid' in bib_data.entries[k].fields.keys():
        title="%s"% \
            (re.sub('[\{\}\\\\]','',bib_data.entries[k].fields['title']))
        #print "Searching for %s in pubmed."%(title)

        result = Entrez.esearch(db=database,term=(title))
        try:
            record=Entrez.read(result)
        except:
            print "Error"
            continue

        print "\n\n{%d} Search for [%s] returned %d records"%(id_num, title,
            int(record['Count']))

        if int(record['Count']) < 1:
            print "Not Found"
        else:
            num=0
            for id in record['IdList']:
                #print id
                title="No Title"
                k_result = Entrez.efetch(db=database, id=id, retmode="xml", \
                    rettype="full", tool="bioPython, Bio.Entrez")
                k_record = Entrez.read(k_result)

                if len(k_record) > 0 and 'MedlineCitation' in k_record[0]:
                    article = k_record[0]['MedlineCitation']['Article']
                    if "ArticleTitle" in article:
                        title = article['ArticleTitle'].encode("utf-8")

                print "[%d] %s %s"%(num,id,title)
                num=num+1

            invar=raw_input("Select the best match or Q for none of these: ")

            if not is_int(invar):
                print "None Found"
            else:
                print "You selected %s: %s"%(invar,record['IdList'][int(invar)])
                bib_data.entries[k].fields['pmid']=record['IdList'][int(invar)]

        w.write_file(bib_data,"out_data_%d.bib"%(id_num))
        id_num=id_num+1

        #result = Entrez.efetch(db=database, id=id, retmode="xml", rettype="full",
		#tool="bioPython, Bio.Entrez")
	#record = Entrez.read(result)
	#if len(record) > 0 and 'MedlineCitation' in record[0]:
		#otherId = record[0]['MedlineCitation']['OtherID']
		#article = record[0]['MedlineCitation']['Article']
		#affiliation = "(blank affiliation)"
		#title = "(blank title)"
		#if "Affiliation" in article:
			#affiliation = article["Affiliation"].encode("utf-8")
		#if "ArticleTitle" in article:
			#title = article['ArticleTitle'].encode("utf-8")
		#print "\t".join([id, title, affiliation]) 
		## print "\t".join([str(i), id, str(openAccess), str(pmc), title.encode("utf-8")])
#
	#else:
		#print "ERROR entry not found for", id
