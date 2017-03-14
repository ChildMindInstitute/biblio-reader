import sys
f = open(sys.argv[1])


print '<?xml version="1.0" encoding="iso-8859-1"?><graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">'

print '<graph edgedefault="undirected" id="">'

for line in f.readlines():
    line = line[:-1]
    line = line.replace('"', "'")
    line = line.replace('<<', "'")
    line = line.replace('>>', "'")
    line = line.replace('<', "")
    line = line.replace('>', "")
    if line[:4] == "node":
        parts = line.split(":")
        print '<node id="' + parts[1] + '" />'
    else:
        parts = line.split(":")
        parts[1] = ":".join(parts[1:])
        terminals = parts[1].split("/")
        if len(terminals) < 2:
            print line
        print '<edge source="' + terminals[0] + '" target="' + terminals[1] + '" />'

print '</graph>'
print '</graphml>'
    
