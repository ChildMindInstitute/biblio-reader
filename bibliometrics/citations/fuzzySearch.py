import os, sys
from difflib import *
import re

junk = None

chars = r"A-Za-z0-9/\-:.,_$%'()[\]<>&; "
shortest_run = 4
regexp = '[%s]{%d,}' % (chars, shortest_run)
pattern = re.compile(regexp)
html_pattern = re.compile(r'<.*?>')
unicode_pattern = re.compile(r'&.*?;')

def process_webarchive(x):
    x = " ".join(pattern.findall(x))
    x = html_pattern.sub('', x)
    x = unicode_pattern.sub('', x)
    ix = x.find("WebResourceResponse_")
    if ix == -1:
        print "PROBLEM!"
        ix = 0
    x = x[:ix]
    return x

def dejunker(x):
    x = x.translate(None, "'\"_.,:").lower().replace(" the "," ").replace("&mdash;", " ").replace("-", " ").replace("  ", " ").replace("\n", "").replace(" ", "")
    return x

titles = []

def getText(file):
    text = open(file, 'r').read()
    if file.endswith("webarchive"):
        text = process_webarchive(text)
    else:
        assert file.endswith("pdf.txt")
    text = dejunker(text)
    return text

def analyze(file):
    print "File :", file
    text = getText(file)
    if len(text) == 0:
        print "ERROR : TEXT LENGTH IS 0"
        return

    m = SequenceMatcher(junk, text, "")
    for title in titles:
        originalTitle = title
        title = dejunker(title)
        ix = text.find(title) 
        ratio = 0.0
        if ix != -1:
            ratio = 1.0
            matchingSubstring = title
            possibleMatch = text[ix:ix+len(title)]
            loc = ix / float(len(text))
        else:
            continue
            m.set_seq2(title)
            result = m.find_longest_match(0, len(text), 0, len(title))
            possibleMatch = text[result.a-result.b:result.a+len(title)-result.b]
            matchingSubstring = text[result.a:result.a+result.size]
            m2 = SequenceMatcher(junk, possibleMatch, title)
            ratio = 0.8 #m2.quick_ratio()
            loc = result.a / float(len(text))
        if ratio > 0.9:
            print "Original title :", originalTitle
            print "Dejunked title :", title
            print "Matching substring :", matchingSubstring
            print "Possible match :", possibleMatch
            print "Ratio :", ratio 
            print "Location :", loc




mendeley = "../data/cmi_library.bib"
lines = open(mendeley, 'r').readlines()
for line in lines:
    line = line[:-1]
    kv = line.split(" = ")
    if kv[0] == "title":
        title = kv[1]
        title = title.replace("{{", "")
        title = title.replace("}},", "")
        titles.append(title)

def walk(numDocs = 1000000):
    root = "../../bin/"
    db = root + "id_to_path.txt"
    i = 0
    ids = []
    filenames = []
    mats = {}
    for line in open(db, 'r').readlines():
        if i == numDocs:
            break
        line = line[:-1]
        if line.startswith("ERROR"):
            continue
        parts = line.split(" : ")
        id = parts[0]
        f = root + parts[1]
        analyze(f)

if __name__ == "__main__":
    walk()
        
            
                          
