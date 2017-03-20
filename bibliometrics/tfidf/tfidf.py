#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
json_to_csv.py

Script to create a csv from the DOT from the Dia diagram.

Authors:
    – Matt Doherty, 2014
    – Jon Clucas, 2017 (jon.clucas@childmind.org)
"""

import math, numpy as np, re, operator, os, textmining
from difflib import *
from collections import defaultdict

rootDir = "../"

def dejunker(x):
    """
    Function to formally standardize words and phrases.
    
    Parameter
    ---------
    x : string
        string to formally standardize
        
    Returns
    -------
    x : string
        the original string `x` formally standardized
    """
    x = x.casefold().replace(" the "," ").replace("-", 
        " ").replace("  ", " ").replace("\n", " ")
    return x

def getText(file):
    """
    Function to get formally standardized text string from a given text file.
    
    Parameter
    ---------
    file : string
        path to file
        
    Returns
    -------
    text : string
        a formally standaradized string of the contents of the text file
        located at `file`.
    """
    # get the text
    text = open(file, 'r').read()
    # standardize the text
    text = dejunker(text)
    # return the standardized text
    return text

def getIDToTitle():
    """
    Function to gather filepaths to scan.
    
    Parameters
    ----------
    None
    
    Returns
    -------
    ret : dictionary
        a dictionary of filepaths
    """
    ret = {}
    refs = os.listdir(os.path.join(rootDir, 'fulltext', 'fulltexts', 'txt'))
    for ref in refs:
        if not ref.startswith('.'):
            ret[ref] = ref
    return ret

def readWordList(filename):
    """
    Function to read a newline delimited textfile.
    
    Parameter
    ---------
    filename : string
        the location of the textfile
        
    Returns
    -------
    ret : set
        set of words, one word per line
    """
    ret = set()
    for line in open(filename,'r').readlines():
        if line.find("#") == -1:
            line = line.rstrip()
            ret.add(line)
    return ret

# define my own tokenizer, so that I can customize it                                                                       
def tokenize(document):
    """
    Function to split strings into words
    
    Parameter
    ---------
    document : string
        a string containing the full text of a document
        
    Returns
    -------
    words : list
        a list of words in the given document
    """
    # make the list case insensitive
    document = document.casefold()
    # regex to find urls
    ur = ('http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-fA-F]'
         '[0-9a-fA-F]))+')
    urls = re.findall(ur, document)
    for i, word in enumerate(urls):
        urls[i] = re.sub('[\)_|]', '', word).strip('.')
    document = re.sub('[^a-z0-9]', ' ', document)
    words = document.strip().split() + urls
    return words

def removeStopwords(words):
    """
    Function to remove unwanted words
    
    Parameter
    ---------
    words : list
        list of words
        
    Returns
    -------
    words : list
        list of words with stopwords removed
    """
    words = [word for word in words if not word.isdigit()]
    words = [word for word in words if word not in textmining.stopwords]
    words = [word for word in words if word not in stopwords]
    return words

# some things have different names, but mean the same thing                                                                     
def standardizeSynonyms(document, dictionary):
    """
    Function to standardize synonyms
    
    Parameters
    ----------
    document : string
        The document to parse synonyms from
       
    dictionary : dictionary
        Dictionary of synonyms with synonyms as keys and canonnical terms as
        values
    """
    words = []
    for w in document:
        if w in dictionary:
            words.append(dictionary[w])
        else:
            words.append(w)
    return words
    
def analyze(file, synDict, stopwords, ngramSetRaw, ngramList):
    """
    Function to get term frequency of a given document
    
    Parameters
    ----------
    file : string
        path to text document
        
    synDict : dictionary
        dictionary of synonyms
        
    stopwords : list
        list of stopwords
        
    ngramSetRaw : set
        set of ngrams
        
    ngramList : list
        list of ngrams
        
    Returns
    -------
    freqs : dictionary
        dictionary of terms : frequencies
    """
    text = getText(file)
    if len(text) == 0:
        print("ERROR : TEXT LENGTH IS 0")
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
    for key in list(counters.keys()):
        if counters[key] >= 15:
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
    """
    Function to calculate TF–IDF for a given term
    
    Parameter
    ---------
    term : string
        a word or phrasal term
    
    Returns
    -------
    tf * idf : float
        term frequency / document frequency
    """
    tf = termFreqs[term]
    idf = math.log(totalDocs / float(docCounters[term]))
    return tf * idf

def calcDf(term):
    """
    Function to calculate DF for a given term
    
    Parameter
    ---------
    term : string
        a word or phrasal term
    
    Returns
    -------
    df : float
        document frequency
    """
    freq = docCounters[term] / float(totalDocs)
    if freq == 0:
        return 0
    else:
        return math.log(freq)

def calcTf(term):
    """
    Function to calculate TF for a given term
    
    Parameter
    ---------
    term : string
        a word or phrasal term
    
    Returns
    -------
    tf : float
        term frequency
    """
    freq = termFreqs[term]
    if freq == 0:
        return 0
    else:
        return 1 + math.log(freq)
                 
if __name__ == "__main__":
    # Synonyms                                                                                       
    synDict = {}
    synSet = readWordList('working/syn_dict.txt')
    for entry in synSet:
        kv = entry.split(',')
        synDict[kv[0]] = kv[1]
    exit
        
    # Stopwords
    stopwords = readWordList('data/stopwords.txt')
    stopwords.update(readWordList('data/extraStopwords.txt'))
    
    # Terms of interest                                                                                      
    searchedTerms = readWordList('data/searchedTerms.txt')

    # Assemble ngrams
    ngramSetRaw = readWordList('working/ngram_list.txt')
    ngramSetRaw.update(searchedTerms)

    ngramList = []
    for term in ngramSetRaw:
        ngramList.append(tuple(term.split()))
    
    i = 0
    idToTitle = getIDToTitle()
    freqsLists = defaultdict(list)
    docCounters = defaultdict(int)
    for id in list(idToTitle.keys()):
        file = os.path.join(rootDir, 'fulltext', 'fulltexts', 'txt', id)
        terms = tokenize(getText(file))
        ngramSetRaw.update(terms)
        for term in terms:
            if term not in ngramList:
                ngramList.append(term)
                
        i += 1
        print(("File :", file))
        docFreqs = analyze(file, synDict, stopwords, ngramSetRaw, ngramList)
        appendDicts(freqsLists, docFreqs)
        mergeDocCounterDicts(docCounters, docFreqs)
    termFreqs = {}
    for key in list(freqsLists.keys()):
        termFreqs[key] = np.median(np.array(freqsLists[key]))
    
    print(("TOTAL DOCS :", i))
    totalDocs = i
    print("MEDIAN TERM FREQS")
    for k in list(termFreqs):
        print((k, ":", termFreqs[k]))
    print("NUM DOCS TERM APPEARED IN")
    for k in docCounters:
        print((k, ":", docCounters[k]))
         
    
    df = {}
    score = {}
    tf = {}
    x = {}
    y = {}
    clazz = {}
    for k in termFreqs:
        if k in docCounters:
            df[k] = calcDf(k)
            tf[k] = calcTf(k)
            score[k] = df[k] * tf[k]
        x[k] = df[k]
        y[k] = tf[k]

    syn = {}
    syn['blood_oxygen_level_dependent'] = "BOLD"
    syn['electroencephalogram'] = "EEG"
    syn['functional_magnetic_resonance_imaging'] = "fMRI"
    syn['prefrontal_cortex'] = "PFC"
    syn['posterior_cingulate'] = "PCC"
    syn['resting_state_network'] = "RSN"
    syn['default_mode_network'] = "DMN"
    syn['positron_emission_tomography'] = "PET"

    labels = {}
    for k in score:
        if k in syn:
            labels[k] = syn[k]
        else:
            if not k in stopwords:
                labels[k] = k
    
    coords = sorted(iter(list(score.items())), key=operator.itemgetter(1))
    coords.reverse()
    
    print(("\t".join(["term", "tf", "df", "score", "x", "y"])))
    with open('data/generated_terms.txt', 'w') as gt:
        for k, v in coords:
            print(("\t".join([labels[k], str(tf[k]), str(df[k]), str(v), str(x[
                  k]), str(y[k])])))
            gt.write(''.join([k, '\n']))
    
          
