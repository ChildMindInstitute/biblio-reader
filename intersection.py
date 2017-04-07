import xml.etree.ElementTree as etree
import pandas as pd
import re

with open('inputs/1000FCP2014.xml', 'r') as f, open('inputs/FCP_DATA.csv') as g:
    data = pd.read_csv(g)
    tree = etree.parse(f)
    root = tree.getroot()
    sente = []
    for name in root.findall(".//{http://www.thirdstreetsoftware" +
                             ".com/SenteXML-1.0}characteristic[@name='articleTitle']"):
        sente.append(str(name.text).lower())

print(sente)
titles = 0
for article in sente:
    if article in data['Title'].apply(lambda x: x.lower()):
        print(article)
        titles += 1
print(titles)

