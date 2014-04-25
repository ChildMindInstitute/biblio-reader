from Bio import Entrez
import textmining
from collections import defaultdict
import re
#from sets import set
	 
Entrez.email='drcc@vt.edu'
database='pubmed'


## alright, first we do an quick analysis of all of the literature on mild traumatic brain injury
neuroQuery = """
("functional neuroimaging" or "positron emission tomography" or PET or magnetoencephalography or MEG or electroencephalography or EEG or "functional magnetic resonance imaging" or fMRI or "diffusion tensor imaging" or DTI or or DSI or "diffusion weighted imaging" or DWI or SWI or "susceptibility weighted imaging" or  "diffusion kurtosis imaging" or "diffusional kurtosis imaging" or DKI or "single photon emission computed tomography" or SPECT or NIRS or "near-infrared spectroscopy" or fNIRS or "functional near-infrared spectroscopy" or "resting state" or "functional connectivity" or "default mode network" )  and ("1995"[Publication Date] : "3000"[Publication Date])
"""

# perform the query the first time to determine
# the number of records
searchResult = Entrez.esearch(db = database, term = neuroQuery)
searchRecord = Entrez.read(searchResult)

# get the number of matches to the query
resultCount = searchRecord['Count']
print "Results:", resultCount

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
