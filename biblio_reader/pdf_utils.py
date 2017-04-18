import pandas as pd
import os
import urllib.request as urllib
import urllib.parse as urlparse
from bs4 import BeautifulSoup as bs

with open('../inputs/FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)

dict_data = dict(zip(data['i'], data['URL']))
valid_data = {key: value for key, value in dict_data.items() if not isinstance(value, float)}
pdfs = []
unlinkable = [5, 7, 52, 59, 95, 116, 304, 341, 345, 359, 385, 386, 478, 491, 495, 500, 624, 653, 666, 702, 756, 794,
              825, 830, 842, 870, 907, 950, 971, 1017, 1020, 1042, 1050, 1093, 1183, 1203, 1219, 1225, 1234, 1235, 1237,
              1287, 1348, 1353, 1355, 1366, 1370, 1372, 1374, 1385, 1386, 1398, 1418, 1431, 1472, 1494, 1546, 1038,
              1202, 1227, 1363, 1373, 1537, 1553, 1556, 127, 282, 364, 413, 490, 524, 758, 796, 844, 961, 1013, 1085,
              1106, 1115, 1146, 1176, 1185, 1211, 1215, 1218, 1221, 1222, 1229, 1233, 1289, 1293, 1305, 1321, 1380,
              1383, 1384, 1393, 1399, 1426, 1449, 1488, 1501]
unlinkable += [key for key, value in valid_data.items() if ('ieeexplore' in value) or ('link.springer' in value) or
               ('liebertpub.com' in value) or ('onlinelibrary.wiley' in value) or ('tandfonline' in value)]
pdfs += [int(pdf.replace('.pdf', '')) for pdf in (os.listdir('../unlink_pdfs') + os.listdir('../linked_pdfs'))] + \
        [key for key, value in valid_data.items() if 'books.google' in value]
no_pdfs = {key: value for key, value in dict_data.items() if key not in pdfs}
unlinkable_data = data[data['i'].isin(unlinkable)]
unlinkable_data.to_csv(path_or_buf='../outputs/unlinkables.csv', index=False)
valid_no_pdfs = {key: value for key, value in dict_data.items() if key not in pdfs and not pd.isnull(value)}

def print_unlinks(unlink_data):
    for row in unlink_data.iterrows():
        row = row[1]
        if 'pdf' not in str(row['URL']):
            print(row['i'])
            print(row['URL'])
            print(row['PMC_LINKS'])


def pdfopener(data):
    for row in data.iterrows():
        url = row[1]['URL']
        if pd.isnull(url):
            continue
        if '.pdf' in url:
            try:
                pdf = urllib.urlopen(url)
                with open('../linked_pdfs/' + str(row[1]['i']) + '.pdf', 'wb') as f:
                    f.write(pdf.read())
            except:
                print(row[1]['i'])
                print(url)
                continue

#def print_unlinkable(data, unlinkable):

def pdffinder(data):
    pdfs_found = 0
    for key in data:
        link = data[key]
        try:
            html = urllib.urlopen(link).read()
        except:
            print(key, link)
            continue
        if isinstance(link, str):
            domain = urlparse.urljoin(link, '/')[:-1]
        else:
            print('No domain name found:', key, type(link))
            continue
        soup = bs(html, 'html.parser')
        for url in soup.find_all('a'):
            if url.get('href') is None:
                continue
            if '.pdf' in url.get('href'):
                pdf = url.get('href')
                if pdf.startswith('/'):
                    pdf = domain + pdf
                try:
                    pdf_file = urllib.urlopen(pdf).read()
                    with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    pdfs_found += 1
                    break
                except Exception as e:
                    print(e, key, pdf)
    print('PDFS FOUND: ' + str(pdfs_found))


def arxiv_open(data):
    for key in data:
        link = data[key]
        try:
            html = urllib.urlopen(link).read()
        except:
            print(key, link)
            continue
        soup = bs(html, 'html.parser')
        for url in soup.find_all('a'):
            if url.get('href') is None:
                continue
            if 'pdf' in url.get('href'):
                pdf = 'https://arxiv.org' + url.get('href') + '.pdf'
                try:
                    pdf_file = urllib.urlopen(pdf).read()
                    with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)

def plos_open(data):
    for key in data:
        link = data[key]
        try:
            html = urllib.urlopen(link).read()
        except:
            print(key, link)
            continue
        soup = bs(html, 'html.parser')
        for url in soup.find_all('a'):
            if url.get('id') is None:
                continue
            if url.get('id') == 'downloadPdf':
                pdf = 'http://journals.plos.org' + url.get('href')
                try:
                    pdf_file = urllib.urlopen(pdf).read()
                    with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)

def liebert_open(data):
    for key in data:
        link = data[key]
        try:
            html = urllib.urlopen(link).read()
        except:
            print(key, link)
            continue
        soup = bs(html, 'html.parser')
        for url in soup.find_all('a'):
            if url.get('href') is None:
                continue
            if 'pdf' in url.get('href'):
                pdf = url.get('href')
                try:
                    pdf_file = urllib.urlopen(pdf).read()
                    with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)

def frontiers_open(data):
    for key in data:
        link = data[key]
        pdf = link.replace('full', 'pdf').replace('abstract', 'pdf')
        if 'pdf' not in pdf:
            pdf = pdf + '/pdf'
        print(pdf)
        try:
            pdf_file = urllib.urlopen(pdf).read()
            with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                f.write(pdf_file)
        except Exception as e:
            print(e, key, pdf)

def citeseer_open(data):
    for key in data:
        link = data[key]
        try:
            pdf_file = urllib.urlopen(link).read()
            with open('../linked_pdfs/' + str(key) + '.pdf', 'wb') as f:
                f.write(pdf_file)
        except Exception as e:
            print(e, key, link)





print(len(no_pdfs))
print(*sorted(valid_no_pdfs.items(), key=(lambda pdf: pdf[1])), sep='\n')
"""
#citeseer_open({key: value for key, value in no_pdfs.items() if 'citeseerx' in value})

print(*[(key, 'http://' + urlparse.urlparse(value).netloc + '.ezproxy.neu.edu' +
         urlparse.urlparse(value).path.replace('full', 'pdf'))
       for key, value in no_pdfs.items() if 'ieee' in value or 'wiley' in value], sep='\n')
"""