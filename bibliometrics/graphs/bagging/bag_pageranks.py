# Usage: bag_pageranks.py original.graphml 1000 0.9

import random, sys, pageranker, operator
random.seed(1)
originalLines = open(sys.argv[1], 'r').readlines()
textFile = 'replicate' + sys.argv[3] + '.txt'
numReplicates = int(sys.argv[2])
prob = float(sys.argv[3])

for i in range(0, numReplicates):
    print "REPLICATE", i
    outFile = open(textFile, 'w')
    for line in originalLines:
        if False: #line.startswith('<edge'):
            if random.uniform(0, 1) < prob:
                outFile.write(line)
        else:
            outFile.write(line)
    outFile.close()
    pageranks = pageranker.getPageRanks(textFile)
    
    sorted_x = sorted(pageranks.iteritems(), key=operator.itemgetter(1))
    for x in sorted_x[-300:]:
        print str(x[0]), ":::", str(x[1])

    
