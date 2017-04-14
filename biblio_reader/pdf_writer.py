import pandas as pd
import re
import os
import urllib.request as urllib
import urllib.parse as urlparse
from bs4 import BeautifulSoup as bs


with open('../inputs/FCP_DATA.csv', 'r', encoding='utf-8') as f:
    data = pd.read_csv(f)


def pdfcollector(data):
    for pdf, authors, title in zip(data['PDF link'], data['Authors'], data['Title']):
        if isinstance(pdf, str):
            pdf = pdf.replace('http://scholar.google.com/', '')
            file = os.path.join('fulltexts', 'pdf', ''.join([re.sub(r'\W', '', authors.split(' & ')[0])[:8], '_',
                                                             re.sub(r'\W', '', title)[:8], '.pdf']))
            if not os.path.exists(file) or os.path.getsize(file) == 0:
                try:
                    pdffile = urllib.request.urlopen(pdf)
                    print(file)
                    with open(file, 'wb') as ofile:
                        ofile.write(pdffile.read())
                except Exception as e:
                    print('\t'.join([str(e), os.path.basename(file), pdf]))
                    continue


def pdffinder(data):
    all_pdfs = []
    for link in data['URL']:
        pdfs = []
        link = str(link).replace('http://scholar.google.com/', '')
        if 'pdf' in link:
            pdfs.append(link)
        else:
            if 'sciencedirect' not in link and 'books.google' not in link and link != 'nan':
                try:
                    html = urllib.urlopen(link)
                except:
                    all_pdfs.append(pdfs)
                    continue
                domain = urlparse.urljoin(link, '/')[:-1]
                soup = bs(html, 'html.parser')
                for url in soup.find_all('a'):
                    if url.get('href') is None:
                        continue
                    if 'pdf' in url.get('href'):
                        pdf = url.get('href')
                        if pdf.startswith('/'):
                            pdf = domain + pdf
                        if pdf not in pdfs:
                            pdfs.append(pdf)
        all_pdfs.append(pdfs)
    data['PDFS'] = all_pdfs

droplist = []
for row in data.iterrows():
    try:
        urllib.urlopen(row[1]['URL'])
        droplist.append(row[0])
    except:
        try:
            urllib.urlopen(row[1]['PMC_LINKS'])
            droplist.append(row[0])
        except:
            continue
no_link_data = data.drop(droplist)


no_link_data.to_csv(path_or_buf='../outputs/unlink_data.csv', index=False)



