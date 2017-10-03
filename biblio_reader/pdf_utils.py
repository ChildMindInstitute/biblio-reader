from bs4 import BeautifulSoup as bs
import os
import pandas as pd
import PyPDF2
import sys
import urllib.parse as urlparse
import urllib.request as urllib
ancestor = os.path.abspath(os.path.join(os.pardir, os.pardir))
if ancestor not in sys.path:
    sys.path.append(ancestor)
from Biblio_Reader import manager as mg


def pdfopener(data, dir):
    """
    Opens all PDF URLs in the dataframe
    """
    for i, row in data.iterrows():
        url = row['URL']
        url = url[26:] if url.startswith(
                'http://scholar.google.com/') else url
        if pd.isnull(url):
            continue
        if '.pdf' in url:
            try:
                pdf = urllib.urlopen(url)
                with open(os.path.join(dir, str(row['i'])) + '.pdf', 'wb') as f:
                    f.write(pdf.read())
            except:
                print('Unable to open:', row['i'])
                print(url)
                continue


def pdffinder(data, dir):
    """
    Tries to find as many PDFs as possible from URLs in the dataframe
    :param data: A zip list of article reference #s and their URLS
    :param dir: Out directory
    """
    data = dict(data)
    pdfs_found = 0
    for key in data:
        link = data[key]
        link = link[26:] if link.startswith(
                'http://scholar.google.com/') else link
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
                    with open(os.path.join(dir, str(key)) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    pdfs_found += 1
                    break
                except Exception as e:
                    print(e, key, pdf)
    print('PDFS FOUND: ' + str(pdfs_found))


def arxiv_open(data, dir):
    """
    Writes all PDFS found on arXiv to the directory
    :param data: A zip list of article reference #s and their URLS
    """
    data = dict(data)
    for key in data:
        link = data[key]
        link = link[26:] if link.startswith(
                'http://scholar.google.com/') else link
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
                    with open(os.path.join(dir, str(key)) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)


def plos_open(data, dir):
    """
    Writes all PDFS found on PLoS website to the dir
    :param data: A zip list of article reference #s and their URLS
    """
    data = dict(data)
    for key in data:
        link = data[key]
        link = link[26:] if link.startswith(
                'http://scholar.google.com/') else link
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
                    with open(os.path.join(dir, str(key)) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)


def liebert_open(data, dir):
    """
    Writes all PDFS found on Liebert website to the directory
    :param data: A zip list of article reference #s and their URLS
    """
    data = dict(data)
    for key in data:
        link = data[key]
        link = link[26:] if link.startswith(
                'http://scholar.google.com/') else link
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
                    with open(os.path.join(dir, str(key)) + '.pdf', 'wb') as f:
                        f.write(pdf_file)
                    break
                except Exception as e:
                    print(e, key, pdf)


def citeseer_open(data, dir):
    """
    Writes all PDFs found on the Citeseerx website to the directory
    :param data: A zip list of article reference #s and their URLS
    """
    data = dict(data)
    for key in data:
        link = data[key]
        link = link[26:] if link.startswith(
                'http://scholar.google.com/') else link
        try:
            pdf_file = urllib.urlopen(link).read()
            with open(os.path.join(dir, str(key)) + '.pdf', 'wb') as f:
                f.write(pdf_file)
        except Exception as e:
            print(e, key, link)


def find_corrupted(pdf_directory):
    """
    Finds all PDFs in the PDF directory that cannot be opened
    :return: List of file names that cannot be opened
    """
    res = []
    for path, dirs, files in os.walk(pdf_directory):
        for file in files:
            full_file = os.path.join(path, file)
            if file == ".DS_Store":
                os.remove(full_file)
            if file.endswith(".pdf"):
                try:
                    PyPDF2.PdfFileReader(full_file)
                except:
                    os.remove(full_file)
                    res.append(int(file.replace('.pdf', '')))
    return res


def find_articles_left(data, dir):
    """
    Finds and prints a list of all articles that do not have PDFs yet
    :param data: A zip list of article reference numbers and their URLS
    :param dir: The PDF directory
    """
    unlinkables = []
    for i, url in data:
        if i not in [file.replace('.txt', '') for file in os.listdir(dir)]:
            print(i, url)
            unlinkables.append((i, url))
    return unlinkables

