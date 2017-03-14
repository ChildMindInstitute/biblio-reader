from nltk.corpus import stopwords
from nltk import tokenize
from collections import *
import operator
from nltk.collocations import *
import nltk
import numpy as np

def doDf():
	counts = Counter()
	i = 0
	allAbstracts = []
	stops = stopwords.words('english')
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			tokens = set(tokenize.wordpunct_tokenize(text))
			tokens = [t for t in tokens if not t in stops]
			for t in tokens:
				counts[t] += 1
	
	sortedCounts = sorted(counts.iteritems(), key=operator.itemgetter(1))
	for elt in sortedCounts[-100:]:
		print elt

def doConditionalTf():
	tfLists = defaultdict(list)
	i = 0
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			if i == 3000:
				break
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			tokens = tokenize.wordpunct_tokenize(text)
			stops = stopwords.words('english')
			tokens = [t for t in tokens if not t in stops]
			tokenSet = set(tokens)
			for t in tokenSet:
				if len(t) > 2:
					tfLists[t].append(text.count(t) / float(len(tokens)))
	tfs = {}
	for t in tfLists.keys():
		tfs[t] = np.median(np.array(tfLists[t]))
			
	sortedTfs = sorted(tfs.iteritems(), key=operator.itemgetter(1))
	for elt in sortedTfs[-1000:]:
		print elt
		
		
def doConditionTfDf():
	counts = Counter()
	i = 0
	allAbstracts = []
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			if i == 3000:
				break
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			allAbstracts.append(text)
	text = " ".join(allAbstracts)
	tokens = tokenize.wordpunct_tokenize(text)
	tokens = [t for t in tokens if not t in stopwords.words('english')]
	for t in tokens:
		counts[t] += 1
	
	sortedCounts = sorted(counts.iteritems(), key=operator.itemgetter(1))
	for elt in sortedCounts[-100:]:
		print elt
		
def doWordFreqs():
	counts = Counter()
	i = 0
	allAbstracts = []
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			if i == 3000:
				break
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			allAbstracts.append(text)
	text = " ".join(allAbstracts)
	tokens = tokenize.wordpunct_tokenize(text)
	tokens = [t for t in tokens if not t in stopwords.words('english')]
	for t in tokens:
		counts[t] += 1
	
	sortedCounts = sorted(counts.iteritems(), key=operator.itemgetter(1))
	for elt in sortedCounts[-100:]:
		print elt

def doBigrams():	
	i = 0
	allAbstracts = []
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			if i == 3000:
				break
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			allAbstracts.append(text)
	text = " ".join(allAbstracts)	
	tokens = tokenize.wordpunct_tokenize(text)
	#tokens = [t for t in tokens if not t in stopwords.words('english')]
	
	bigram_measures = nltk.collocations.BigramAssocMeasures()
	trigram_measures = nltk.collocations.TrigramAssocMeasures()
	
	# change this to read in your data
	finder = BigramCollocationFinder.from_words(tokens)
	
	# only bigrams that appear 3+ times
	finder.apply_freq_filter(200) 
	
	# return the 10 n-grams with the highest PMI
	ret = finder.nbest(bigram_measures.pmi, 20)  
	
	for elt in ret:
		print elt
		
def doTrigrams():	
	i = 0
	allAbstracts = []
	for line in open('../data/new_pmids.bib').readlines():	
		line = line.rstrip()
		if line.startswith("abstract"):
			i += 1
			text = line.split(" = ")[1]
			text = text[1:-2] # remove leading " and trailing ",
			text = text.lower()
			allAbstracts.append(text)
	text = " ".join(allAbstracts)	
	tokens = tokenize.wordpunct_tokenize(text)
	#tokens = [t for t in tokens if not t in stopwords.words('english')]
	
	bigram_measures = nltk.collocations.BigramAssocMeasures()
	trigram_measures = nltk.collocations.TrigramAssocMeasures()
	
	# change this to read in your data
	finder = TrigramCollocationFinder.from_words(tokens)
	
	# only bigrams that appear 3+ times
	finder.apply_freq_filter(100) 
	
	# return the 10 n-grams with the highest PMI
	ret = finder.nbest(trigram_measures.pmi, 20)  
	
	for elt in ret:
		print elt
		
		
if __name__ == "__main__":
	doDf()