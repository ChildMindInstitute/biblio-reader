import sys, unicodedata
from unidecode import unidecode
import HTMLParser
h = HTMLParser.HTMLParser()

f = open(sys.argv[1], 'r')
ofile = open(sys.argv[2], 'w')
lines = f.readlines()
for line in lines:
    #line = line.decode('utf-8')
    line = line.decode('utf-8')
    #line = unicodedata.normalize('NFKD', line)
    line = h.unescape(line)
    line = unidecode(line)
    line = line.replace("&", "")
    line = line.replace("#", "")
    line = line.encode('ascii')
    ofile.write(line)

ofile.close()
