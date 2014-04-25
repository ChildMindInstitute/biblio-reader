# -*- coding: utf-8 -*-
import os, sys
import re
from DbAPI import DbAPI 
import os.path
import unidecode, unicodedata

rootDir = "/home/files/txt/"

def dejunker(x):
	if not isinstance(x, unicode):
		x = unicode(x, "utf-8")
		
	x = x.replace("The phase shift index for marking functional asynchrony in Alzheimer's disease patients using fMRI.", "The phase shift index for marking functional asynchrony in Alzheimer's patients by fMRI.")
	x = x.replace("Residual functional connectivity in the split-brain revealed with resting-state functional MRI.", "Residual functional connectivity in the split-brain revealed with resting-state fMRI.")
	x = x.replace("Neurophysiological responses to traumatic reminders in the acute aftermath of serious motor vehicle collisions using [15O]-H2O positron emission tomography.", "Neurophysiological responses to traumatic reminders in the acute aftermath of serious motor vehicle collisions using [15O]-H2O PET.")
	x = x.replace("Modulation of a human memory circuit by subsyndromal depression in late life: a functional magnetic resonance imaging study.", "Modulation of a human memory circuit by subsyndromal depression in late life: an fMRI study.")
	x = x.replace("The ", "").lower().replace("\n", " ").replace(" the ", "").replace(" of ", "").replace(" ", "").replace(u"ﬁ", "fi").replace("$", "").replace(",", "").replace(".", "").replace(u"—", "-").replace(u"–", "-").replace("-", "").replace(":", "").replace(u"’", "'").replace(u"‘", "'").replace("''", "'").replace(u"“", "'").replace(u"”", "'").replace('"', "'").replace("modelling", "modeling")
	#x = unidecode.unidecode(x)
	#if not isinstance(x, unicode):
	#	x = unicode(x, "utf-8", errors="replace")
	x = unicodedata.normalize('NFKD', x)	
	x = x.replace("decreaseddefaultmodeneuralmodulationwithageinschizophrenia", "decreaseddefaultmodeneuralnetworkmodulationwithageinschizophrenia")
	x = x.replace("selffacerecognitionactivatesafrontoparietalbmirrorqnetworkinrighthemisphereaneventrelatedfmristudy", "selffacerecognitionactivatesafrontoparietal'mirror'networkinrighthemisphereaneventrelatedfmristudy")
	x = x.replace("actionrelatedpropertiesshapeobjectrepresentationsinventralstream", "actionrelatedpropertiesobjectsshapeobjectrepresentationsinventralstream")
	x = x.replace("functionalmriinadhdasystematicliteraturereview", "functionalmagneticresonanceimaginginattentiondeficithyperactivitydisorder(adhd)asystematicliteraturereview")
	x = x.replace("dynamicadjustmentsinprefrontalhippocampalandinferiortemporalinteractionswithincreasingvisualworkingmemoryload", "dynamicadjustmentsinfrontalhippocampalandinferiortemporalinteractionswithincreasingvisualworkingmemoryload")
	return x

def getText(file):
	text = open(file, 'r').read()
	text = dejunker(text)
	return text

def analyze(file, trueTitle, titles):
	print "File :", file
	
	text = getText(file)
	if len(text) < 1000:
		print "ERROR : TEXT NOT LONG ENOUGH"
		return
	
	for originalTitle in titles:
		potTitle = dejunker(originalTitle)
		ix = text.find(potTitle)
		
		selfTitle = False 
		if potTitle == dejunker(trueTitle):
			selfTitle = True
			if ix == -1:
				#print "Dejunked true title:"
				#print dejunker(trueTitle)
				#print "Pot title:"
				#print potTitle
				#print text
				#print type(text)
				#print type(potTitle)
				print "ERROR : Couldn't find title '%s' in source '%s'" % (originalTitle.encode('ascii', 'ignore'), file)
			
		if ix != -1 or selfTitle:
			ratio = 1.0
			matchingSubstring = potTitle
			possibleMatch = text[ix:ix+len(potTitle)]
			loc = ix / float(len(text))
			print "Original title :", originalTitle.encode('ascii', 'ignore')
			print "Dejunked title :", potTitle.encode('ascii', 'ignore')
			print "Matching substring :", matchingSubstring.encode('ascii', 'ignore')
			print "Possible match :", possibleMatch.encode('ascii', 'ignore')
			print "Self title :", selfTitle
			print "Ratio :", ratio 
			print "Location :", loc
	print "Finished :", file

def getIDToTitle():	
	db = DbAPI()
	ret = {}
	refs = db.execute("select * from refs")
	for ref in refs:
		ret[ref['pmid']] = ref['title']
	return ret
		
if __name__ == "__main__":
	idToTitle = getIDToTitle()
	ids = idToTitle.keys()
	ids.sort()
	for id in ids:
		f = rootDir + str(id) + ".pdf.txt"
		if os.path.exists(f):
			analyze(f, idToTitle[id], idToTitle.values())
