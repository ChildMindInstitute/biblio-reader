import xml.etree.ElementTree as etree
import pandas as pd

with open('inputs/1000FCP2014.xml', 'r') as f, open('inputs/FCP_DATA.csv') as g:
    data = pd.read_csv(g)
    tree = etree.parse(f)
    root = tree.getroot()
    sente = []
    for name in root.findall(".//{http://www.thirdstreetsoftware" +
                             ".com/SenteXML-1.0}characteristic[@name='articleTitle']"):
        sente.append(str(name.text).lower())

titles = 0
for row in data.iterrows():
    if row[1]['Title'].lower() in sente and row[1]['Year'] > 2014:
        print(str(row[1]['Title']) + "   " + str(int(row[1]['Year'])) + "    " + str(row[0]))
        titles += 1
print(titles)

