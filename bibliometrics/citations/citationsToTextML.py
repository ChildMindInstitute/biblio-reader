# Usage : self.py searchCitationsOutput.txt > output.textml
import sys

f = open(sys.argv[1], 'r')
lines = f.readlines()
nodes = set()
edgesList = []

for line in lines:
	line = line[:-1]
	if ":" in line:
		parts = line.split(" : ")
		if parts[0] == "Finished":
			if startPoint == None:
				assert False, "PROBLEM WITH %s" % currentFile
				break
			for endPoint in endPoints:
				edgesList.append((startPoint, endPoint))			
		elif parts[0] == "File":
			currentFile = parts[1]	
			endPoints = []
			startPoint = None
			currentTitle = None
			trueTitle = None
			currentRatio = None			
			selfTitle = None			
		elif parts[0] == "Original title":
			nodes.add(parts[1])
			currentTitle = parts[1]
		elif parts[0] == "Ratio":
			currentRatio = float(parts[1])
		elif parts[0] == "Self title":
			if parts[1] == "True":
				selfTitle = True
			else:
				selfTitle = False
		elif parts[0] == "Location":
			location = float(parts[1])
			if currentRatio > 0.9:
				if selfTitle:
					startPoint = currentTitle
				else:
					endPoints.append(currentTitle)
			selfTitle = None
			currentRatio = None
			currentTitle = None

for node in nodes:
	print 'node:' + node
for edge in edgesList:
	print 'edge:' + edge[0] + "/" + edge[1]



