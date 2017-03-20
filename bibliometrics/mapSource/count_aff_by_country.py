import sys
from collections import Counter
import operator

countries = []
for line in open('countries.txt', 'r').readlines():
	line = line.rstrip()
	line = line.replace('"', '').lower()
	countries.append(line)
countrySet = set(countries)
	
def getCountryCounts():
	counts = Counter()
	for line in open('../data/cmi_library.bib', 'r').readlines():
		line = line[:-1]
		if "address = " in line:
			kv = line.split(" = ")
			country = kv[1].split("{")[1].split("}")[0].lower()
			if country == "england":
				country = "united kingdom"
			if "new york" in country:
				country = "united states"
			if "cambridge" in country:
				country = "united kingdom"
			if "korea" in country:
				country = "korea, republic of"
			if country not in countrySet:
				for pot in countries:
					if pot in country or country in pot:
						# print "Matched", country, "with", pot
						country = pot
						break
			if country not in countrySet:
				print "Couldn't figure out", country
			counts[country] += 1
			
	return counts
	
	
def getAffiliationCounts(): 
	corrections = [("jersey", "united states"), ("georgia", "united states"), ("belgique", "belgium"), ("hunan", "china"), ("chongqing", "china"), ("stockholm", "sweden"), ("frankfurt", "germany"), ("dallas", "united states"), ("nashville", "united states"), ("michigan", "united states"), ("baltimore", "united states"), ("philadelphia", "united states"), ("new haven", "united states"), ("iran", "iran (islamic republic of)"), ("washington", "united states"), ("tokyo", "japan"), ("kyoto", "japan"), ("deutschland", "germany"), ("korea", "korea, republic of"), (", uk", "united kingdom"), ("england", "united kingdom"), ("new york", "united states"), ("usa", "united states"), ("cambridge", "united kingdom")]
	removes = set(["jordan", "jersey", "georgia"])
	counts = Counter()
	lines = open('affiliations.txt', 'r').readlines()
	for line in lines:
		line = line.rstrip().lower()
		if "found 1614 results" in line or "error" in line or "blank affiliation" in line:
			continue
		line = line.split("\t")[2]
		something = False
		for pot in countries:
			if pot in line and pot not in removes:
				counts[pot] += 1
				something = True
		for correction in corrections:
			if correction[0] in line and correction[1] not in line:
				counts[correction[1]] += 1
				something = True
		# if not something:
			# print "No country for", line
	return counts
				
			
counts = getAffiliationCounts()

print "count"
covered = set()
for country in countries:
	covered.add(country)
	print counts[country]
	
for country in counts.keys():
	if not country in covered:
		print "I have a count for", country, "but nowhere to put it"
  
sorted_x = sorted(counts.iteritems(), key=operator.itemgetter(1))
total = 0
for x in sorted_x:
	# print x
	total += x[1]
# print "total", total
