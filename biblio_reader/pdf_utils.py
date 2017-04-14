import pandas as pd
import os
import urllib.request as urllib

with open('../outputs/unlink_data.csv', 'r') as f, open('../inputs/FCP_DATA.csv', 'r') as g, \
        open('../outputs/FCP_DATA_PDFS.csv', 'r') as h:
    unlink_data = pd.read_csv(f)
    data = pd.read_csv(g)
    pdf_data = pd.read_csv(h)


"""
for row in unlink_data.iterrows():
    row = row[1]
    if 'pdf' not in str(row['URL']):
        print(row['Unnamed: 0'])
        print(row['URL'])
        print(row['PMC_LINKS'])

"""

unlinkable = [95, 116, 304, 341, 345, 359, 386, 478, 491, 495, 500, 624, 653, 666, 702, 756, 794, 825, 830, 842, 870,
              907, 950, 971, 1017, 1020, 1042, 1050, 1093, 1183, 1203, 1219, 1225, 1234, 1235, 1237, 1287, 1348, 1353,
              1355, 1366, 1370, 1372, 1374, 1385, 1386, 1398, 1418, 1431, 1472, 1494, 1546, 1038, 1202, 1227, 1363]

pdfs = []

for pdf in os.listdir('../unlink_pdfs'):
    pdfs.append(pdf.replace('.pdf', ''))

for pdf in os.listdir('../linked_pdfs'):
    pdfs.append(pdf.replace('.pdf', ''))

print(len(pdfs))


def pdfopener(data):
    for row in data.iterrows():
        url = row[1]['URL']
        if pd.isnull(url):
            continue
        if '.pdf' in url:
            try:
                pdf = urllib.urlopen(url)
                with open('../linked_pdfs/' + str(row[1]['Unnamed: 0']) + '.pdf', 'wb') as f:
                    f.write(pdf.read())
            except:
                print(row[1]['Unnamed: 0'])
                print(url)
                continue
