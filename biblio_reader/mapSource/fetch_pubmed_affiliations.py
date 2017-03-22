from Bio import Entrez
from  pybtex.database.input import bibtex
import textmining, pprint, csv
from collections import defaultdict
import re, sys
import pprint

Entrez.email = 'drcc@vt.edu'
database = 'pubmed'

#lines = open('../data/cmi_library.bib', 'r').readlines()
#for line in lines:
	#line = line.rstrip()
	#if line.startswith("pmid = "):
		#pmid = line.split(" = ")[1][1:-2]
		#pmids.append(pmid)

parser=bibtex.Parser()
bib_data=parser.parse_file('../data/cmi_library.bib')

pmids = []
for k in bib_data.entries.keys():
    if not 'pmid' in bib_data.entries[k].fields.keys():
        title="%s"% \
            (re.sub('[\{\}\\\\]','',bib_data.entries[k].fields['title']))
        print "Key %s is missing: %s"%(k,title)
    else:
        if bib_data.entries[k].fields['pmid'] not in pmids:
             pmids.append(bib_data.entries[k].fields['pmid'])
        else:
            print "%s is a duplicate"%(bib_data.entries[k].fields['pmid'])

print "Found", len(pmids), "results"

i = 0
total = len(pmids)

pp = pprint.PrettyPrinter(indent=4)
for id in pmids:
	i += 1

	if i < -1:
		continue
	id = id.rstrip()
	# get the record

	result = Entrez.esearch(db=database, term=(id))
	try:
		record = Entrez.read(result)
	except:
		print "ERROR"
		continue

	# print "OpenAccess:", openAccess
	
	result = Entrez.efetch(db=database, id=id, retmode="xml", rettype="full",
		tool="bioPython, Bio.Entrez")
	record = Entrez.read(result)
	if len(record) > 0 and 'MedlineCitation' in record[0]:
		otherId = record[0]['MedlineCitation']['OtherID']
		article = record[0]['MedlineCitation']['Article']
		affiliation = "(blank affiliation)"
		title = "(blank title)"
		if "Affiliation" in article:
			affiliation = article["Affiliation"].encode("utf-8")
		if "ArticleTitle" in article:
			title = article['ArticleTitle'].encode("utf-8")
                print "Affiliation: %s"%(affiliation)
                print "Title: %s"%(title)
		print "\t".encode("utf-8").join([id.encode("utf-8"), title, affiliation]) 
		#print "\t".join([str(i), id, str(openAccess), str(pmc), title.encode("utf-8")])

	else:
		print "ERROR entry not found for", id
