import sys, operator, numpy

titleToRanks = {}

for line in open(sys.argv[1], 'r').readlines():
    line = line[:-1]
    if ":::" in line:
        parts = line.split(" ::: ")
        title = parts[0].replace("\t", " ")
        rank = float(parts[1])
        if not title in titleToRanks:
            titleToRanks[title] = []
        titleToRanks[title].append(rank)

means = {}
stddevs = {}

for title, ranks in titleToRanks.items():
    means[title] = numpy.mean(ranks)
    stddevs[title] = numpy.std(ranks)

allSum = 0
for title in means.keys():
    allSum += float(means[title])

sorted_x = sorted(means.iteritems(), key=operator.itemgetter(1))

print "title\tmean\tstddev"
topSum = 0
for x in sorted_x[-10:]:
    topSum += x[1]
    print "\t".join([str(x[0]), str(x[1]), str(stddevs[x[0]])])

print topSum, allSum, (topSum / allSum)
