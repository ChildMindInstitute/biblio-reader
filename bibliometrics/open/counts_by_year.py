from Bio import Entrez
import textmining
from collections import defaultdict
import re, sys
#from sets import set
	 
pubmedFile = open('../data/open_rs_pubmed.txt', 'r')
thisEntry = None
pubmedEntries = []
for line in pubmedFile:
        line = line[:-1]
	if not ": " in line:
		continue
	if line.startswith("Counter"):
                if thisEntry != None:
                        pubmedEntries.append(thisEntry)
		thisEntry = {}
        else:
                parts = line.split(": ")
                thisEntry[parts[0]] = parts[1]

rsFreeByYear = defaultdict(int)
for entry in pubmedEntries:
	if entry['OpenAccess'] == "True":
		rsFreeByYear[entry['Year']] += 1

rsByYear = defaultdict(int)
for line in open('../data/cmi_library.bib').readlines():
	if line[0] == "@":
		name = line.split("{")[1][:-1]
	else:
		kv = line.split(" = ")
		if kv[0] == "year":
			value = kv[1][1:-2].title()
			rsByYear[value] += 1

Entrez.email='drcc@vt.edu'
database='pubmed'

## alright, first we do an quick analysis of all of the literature on mild traumatic brain injury
neuroQuery = """
("fMRI" or "functional magnetic resonance imaging") and ("YEAR_HERE"[Publication Date])
"""

print "\t".join(["year", "neuroTotal", "neuroFree", "rsTotal", "rsFree"])
for year in range(1995, 2013):
	year = str(year)
	query = neuroQuery.replace("YEAR_HERE", year)
	searchResult = Entrez.esearch(db = database, term = query)
	searchRecord = Entrez.read(searchResult)

	# get the number of matches to the query
	neuroTotal = searchRecord['Count']
	
	query += " AND free full text[sb]"
	searchResult = Entrez.esearch(db = database, term = query)
	searchRecord = Entrez.read(searchResult)
	neuroFree = searchRecord['Count']

	print "\t".join([year, neuroTotal, neuroFree, str(rsByYear[year]), str(rsFreeByYear[year])])





sys.exit(0)

# perform the query again to get all of the PMIDs that match
searchResult = Entrez.esearch(db = database, term = neuroQuery, retmax = resultCount)
searchRecord = Entrez.read(searchResult)


for id in searchRecord['IdList']:
	print id
	continue
	#get the record
	id_dl_handle=Entrez.efetch(db=database,id=id,retmode="xml",rettype="full",
		tool="bioPython, Bio.Entrez")
	id_dl_record=Entrez.read(id_dl_handle)

	if 'MedlineCitation' in id_dl_record[0]:
		# add in the publication year to the pub year hist
		try:
			pub_year=id_dl_record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
		except:
			pub_year=id_dl_record[0]['MedlineCitation']['DateCreated']['Year']
	
		query2_year_hist[pub_year]+=1
	
		# add the journal name into the journal name hist
		journal_abbrev=id_dl_record[0]['MedlineCitation']['Article']['Journal']['Title']
		query2_journal_hist[journal_abbrev]+=1

		# add in any meshterms that might exist
		if "MeshHeadingList" in id_dl_record[0]['MedlineCitation']:
			mesh_terms_list=id_dl_record[0]['MedlineCitation']['MeshHeadingList']
			for mt in mesh_terms_list:
				mesh_term=mt['DescriptorName'].title()
				query2_mesh_hist[mesh_term]+=1

		# add in funding agencies into hist
		if 'GrantList' in id_dl_record[0]['MedlineCitation']['Article']:
			grantlist=id_dl_record[0]['MedlineCitation']['Article']['GrantList']
			for grant in grantlist:
				if 'Agency' in grant:
					agency=grant['Agency']
					query2_agency_hist[agency]+=1

		# add in title keywords
		if 'ArticleTitle' in id_dl_record[0]['MedlineCitation']['Article']:
			title=id_dl_record[0]['MedlineCitation']['Article']['ArticleTitle']
			#if title.find('cerebral') != -1:
				#print title
			#title_words=textmining.simple_tokenize_remove_stopwords(title)
			title_words=cc_tokenize_remove_stopwords(title)
			title_words=textmining.collapse_ngrams(title_words,ngram_list)	
			title_words=cc_standardize_synonyms(title_words,syn_dict)
			#title_words=textmining.stem(title_words)
			title_words=set(title_words) # make unique
			if "magnetic_resonance" in title_words:
				print title
			for w in title_words:
				query2_title_keywords_hist[w] += 1
			title_keywords.add_doc(id_dl_record[0]['MedlineCitation']['Article']['ArticleTitle'])

		# add in abstract keywords
		if 'Abstract' in id_dl_record[0]['MedlineCitation']['Article']:
			abstract=id_dl_record[0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
			#if abstract.find('cerebral') != -1:
				#print abstract
			#abstract_words=textmining.simple_tokenize_remove_stopwords(abstract)
			abstract_words=cc_tokenize_remove_stopwords(abstract)
			abstract_words=textmining.collapse_ngrams(abstract_words,ngram_list)	
			abstract_words=cc_standardize_synonyms(abstract_words,syn_dict)
			#abstract_words=textmining.stem(abstract_words)
			if "magnetic_resonance" in abstract_words:
				print abstract
			abstract_words=set(abstract_words)
			for w in abstract_words:
				query2_abstract_keywords_hist[w] += 1
			abstract_keywords.add_doc(id_dl_record[0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0])

		# combine keywords for title and abstract
		keywords=title_words.union(abstract_words)
		for w in keywords:	
			query2_keywords_hist[w] += 1

	else:
		print "Medline not found", id_dl_record[0].keys()

title_keywords.write_csv('query2_title_keywords.csv',cutoff=1);
abstract_keywords.write_csv('query2_abstract_keywords.csv',cutoff=1);

ofile=open('query2_year_hist.csv','w')
for items in query2_year_hist.items():
	ofile.write("%s,%s\n" % items)
ofile.close

keywords=query2_keywords_hist.keys()
keywords.sort()
ofile=open('query2_keywords_hist.csv','w')
for kw in keywords:
	if (int(query2_keywords_hist[kw])>1):
		ofile.write("%s,%s\n" % (kw,query2_keywords_hist[kw]))
ofile.close

ofile=open('query2_imaging_words_hist.csv','w')
for kw in keywords:
	if kw in imaging_words:
		ofile.write("%s,%s\n" % (kw,query2_keywords_hist[kw]))
ofile.close


ofile=open('query2_title_words_hist.csv','w')
for items in query2_title_keywords_hist.items():
	if(int(items[1])>1):
		ofile.write("%s,%s\n" % items)
ofile.close

ofile=open('query2_abstract_words_hist.csv','w')
for items in query2_abstract_keywords_hist.items():
	if(int(items[1])>1):
		ofile.write("%s,%s\n" % items)
ofile.close

ofile=open('query2_mesh_hist.csv','w')
for items in query2_mesh_hist.items():
	ofile.write("%s,%s\n" % items)
ofile.close

ofile=open('query2_agency_hist.csv','w')
for items in query2_agency_hist.items():
	ofile.write(u'%s,%s\n' % items)
ofile.close

ofile=open('query2_journal_hist.csv','w')
for items in query2_journal_hist.items():
	try:
		ofile.write(u'%s,%s\n' % items)
	except:
		continue
ofile.close
