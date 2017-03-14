import textmining
import os, sys, operator
from difflib import *
from collections import defaultdict
import numpy as np
import DbAPI

rootDir = "/home/files/txt/"

def dejunker(x):
	x = x.translate(None, "'\"_.,:").lower().replace(" the "," ").replace("-", " ").replace("  ", " ").replace("\n", " ").replace(" ", "")
	return x

def getText(file):
	text = open(file, 'r').read()
	text = dejunker(text)
	return text

def getIDToTitle():	
	db = DbAPI()
	ret = {}
	refs = db.execeute("select * from refs")
	for ref in refs:
		ret[ref['pmid']] = ref['title']
	return ret

def readWordList(filename):
	ret = set()
	for line in open(filename,'r').readlines():
		if line.find("#") == -1:
			line = line.rstrip()
			ret.add(line)
	return ret

# define my own tokenizer, so that I can customize it																	   
def tokenize(document):
	document = document.lower()
	document = re.sub('[^a-z0-9]', ' ', document)
	words = document.strip().split()
	return words

def removeStopwords(words):
	words = [word for word in words if word not in textmining.stopwords]
	words = [word for word in words if word not in stopwords]
	words = [word for word in words if word not in extraStopwords]
	return words

# some things have different names, but mean the same thing																     
def standardizeSynonyms(document, dictionary):
	words = []
	for w in document:
		if w in dictionary:
			words.append(dictionary[w])
		else:
			words.append(w)
	return words
	
def analyze(file, synDict, stopwords, ngramSetRaw, ngramList):
	text = getText(file)
	if len(text) == 0:
		print "ERROR : TEXT LENGTH IS 0"
		return {}

	tokens = tokenize(text)
	totalWords = len(tokens)
	words = removeStopwords(tokens)
	words = textmining.collapse_ngrams(words, ngramList)
	words = standardizeSynonyms(words, synDict)
	counters = defaultdict(int)
	for word in words:
		if word in ngramSetRaw:
			counters[word] += 1
	freqs = {}
	for key in counters.keys():
		freqs[key] = counters[key] / float(totalWords)
	return freqs


def appendDicts(a, b):
	if len(b) > 0:
		for k in b:
			a[k].append(b[k])

def mergeDocCounterDicts(a, b):
	if len(b) > 0:
		for k in b:
			a[k] += 1

def calcTfidf(term):
    tf = termFreqs[term]
    idf = math.log(totalDocs / float(docCounts[term]))
    return tf * idf

def calcDf(term):
    freq = docCounts[term] / float(totalDocs)
    if freq == 0:
	    return 0
    else:
	    return math.log(freq)

def calcTf(term):
    freq = termFreqs[term]
    if freq == 0:
	    return 0
    else:
	    return 1 + math.log(freq)

def getClazz(k, imagingTerms, cogAtlasTerms, pubBrainTerms):
	if k in imagingTerms:
		clazz[k] = "Imaging"
	elif k in cogAtlasTerms:
		clazz[k] = "Cognitive Atlas"
	elif k in pubBrainTerms:
		clazz[k] = "PubBrain"
	else:
		clazz[k] = "Librarian Tags"
					
if __name__ == "__main__":
	# Synonyms																				       
	synDict = {}
	synSet = readWordList('data/synDict.txt')
	for entry in synSet:
		kv = entry.split(',')
		synDict[kv[0]] = kv[1]
		
	# Stopwords
	stopwords = readWordList('data/stopwords.txt')
	stopwords.update(readWordList('data/extraStopwords.txt'))
	
	# Terms of interest																					  
	imagingTerms = readWordList('data/imagingTerms.txt')
	cogAtlasTerms = readWordList('data/cogAtlas.csv')
	pubBrainTermsRaw = readWordList('data/pubBrainLexicon.txt')
	pubBrainTerms = set()
	for term in pubBrainTermsRaw:
		term = term.replace("*","")
		pubBrainTerms.add(term) 
	
	# Assemble ngrams
	ngramSetRaw = readWordList('data/ngramList.txt')
	ngramSetRaw.update(imagingTerms)
	ngramSetRaw.update(cogAtlasTerms)
	ngramSetRaw.update(pubBrainTerms)
	ngramList = []
	for term in ngramSetRaw:
		ngramList.append(tuple(term.split()))
	
	i = 0
	idToTitle = getIDToTitle()
	freqsLists = defaultdict(list)
	docCounters = defaultdict(int)
	for id in idToTitle.keys():
		file = rootDir + id + ".txt"
		i += 1
		print "File :", file
		docFreqs = analyze(file, synDict, stopwords, ngramSetRaw, ngramList)
		appendDicts(freqsLists, docFreqs)
		mergeDocCounterDicts(docCounters, docFreqs)
	
	termFreqs = {}
	for key in freqsLists.keys():
		termFreqs[key] = np.median(np.array(freqsLists[key]))
	
	print "TOTAL DOCS :", i
	print "MEDIAN TERM FREQS"
	for k in termFreqs:
	    print k, ":", termFreqs[k]
	print "NUM DOCS TERM APPEARED IN"
	for k in docCounters:
	    print k, ":", docCounters[k]
		 
	
	df = {}
	score = {}
	tf = {}
	x = {}
	y = {}
	clazz = {}
	for k in termFreqs:
	    if k in docCounts:
	        df[k] = calcDf(k)
	        tf[k] = calcTf(k)
	        score[k] = df[k] * tf[k]
		x[k] = df[k]
		y[k] = tf[k]
		clazz[k] = getClazz(k, imagingTerms, cogAtlasTerms, pubBrainTerms)
	
	tweak = {}
	tweak['insula'] = (0, 0.09)
	tweak['precuneus'] = (0, 0.09)
	tweak['thalamus'] = (0, 0.09)
	tweak['encoding'] = (0, 0.09)
	tweak['cognition'] = (0.025, 0.09)
	tweak['positron_emission_tomography'] = (-0.01, 0.09)
	tweak['anterior_cingulate'] = (0.075, 0)
	tweak['activation'] = (-0.02, -0.09)
		

	
	for k in y.keys():
		if not k in tweak:
			y[k] -= 0.09
	
	for k in tweak:
	 	x[k] += tweak[k][0]
	 	y[k] += tweak[k][1]


	syn = {}
	syn['blood_oxygen_level_dependent'] = "BOLD"
	syn['electroencephalogram'] = "EEG"
	syn['functional_magnetic_resonance_imaging'] = "fMRI"
	syn['prefrontal_cortex'] = "PFC"
	syn['posterior_cingulate'] = "PCC"
	syn['resting_state_network'] = "RSN"
	syn['default_mode_network'] = "DMN"
	syn['positron_emission_tomography'] = "PET"

	remove = {}
	# ?remove['hippocampus'] = True
	# remove['thalamus'] = True
	# remove['language'] = True
	remove['resting_state'] = True
	remove['functional_magnetic_resonance_imaging'] = True
	remove['functional_connectivity'] = True
	remove['connectivity'] = True
	remove['search'] = True
	remove['blood_oxygen_level_dependent'] = True
		
	labels = {}
	for k in score:
		if k in syn:
			labels[k] = syn[k]
		else:
			labels[k] = k
	
	coords = sorted(score.iteritems(), key=operator.itemgetter(1))
	coords.reverse()
	
	print "\t".join(["term", "tf", "df", "score", "x", "y", "Source"])
	for k, v in coords:
		if not k in remove:
			print "\t".join([labels[k], str(tf[k]), str(df[k]), str(v), str(x[k]), str(y[k]), clazz[k]])
	 
	
	      
