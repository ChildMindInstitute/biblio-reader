import xml.etree.ElementTree as etree
import pandas as pd

with open('inputs/1000FCP2014.xml', 'r') as f, open('inputs/FCP_DATA.csv') as g:
    data = pd.read_csv(g)
    tree = etree.parse(f)

root = tree.getroot()
sente = []
for name in root.findall(".//{http://www.thirdstreetsoftware.com/SenteXML-1.0}characteristic[@name='articleTitle']"):
    sente.append(str(name.text).lower())
"""
titles = 0
for article in sente:
    if article in data['Title'].apply(lambda x: x.lower()):
        print(article)
        titles += 1
print(titles)
"""

authors = []
for author in root.findall(".//{http://www.thirdstreetsoftware.com/SenteXML-1.0}authors"):
    authors.append(' '.join(author.itertext()).strip().split('\n'))

sentezip = zip(sente, authors)

def make_dataframe(zip):
    list = []
    for sente, authors in zip:
        list.append([sente, authors])
    return pd.DataFrame(data=list, columns=['sente', 'authors'])

make_dataframe(sentezip).to_csv(path_or_buf='outputs/intersections.csv')