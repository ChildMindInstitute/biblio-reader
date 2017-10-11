from bs4 import BeautifulSoup as bs
import os
import pandas as pd
import PyPDF2
from pdfrw import PdfReader
import sys
import urllib.parse as urlparse
import urllib.request as urllib
ancestor = os.path.abspath(os.path.join(os.pardir, os.pardir))
if ancestor not in sys.path:
    sys.path.append(ancestor)
from Biblio_Reader import manager as mg
from numpy import nan


def check_url(s):
    """
    Function to check for existance and 'http://scholar.google.com/' prefix
    
    Parameters
    ----------
    s : string
    
    Returns
    -------
    s : string or None
    """
    if ((s in {"", nan, "nan", None}) or
        (len(s) <= 3)):
        return(None)
    elif s.startswith(
        'http://scholar.google.com/'):
        return(s[26:])
    else:
        return(s)
        

def pdf_anchor_finder(url):
    """
    Function to find pdf links on a given page
    
    Parameters
    ----------
    url : string
    
    Returns
    -------
    pdf_urls : list of strings
    """
    try:
        html = urllib.urlopen(url).read()
    except:
        return([])
    if isinstance(url, str):
        domain = urlparse.urljoin(url, '/')[:-1]
    else:
        return([])
    soup = bs(html, 'html.parser')
    pdf_urls = [url] if 'pdf' in url else []
    for url_a in soup.find_all('a'):
        if url_a.get('href') is None:
            return([])
        if 'pdf' in url_a.get('href'):
            pdf = url_a.get('href')
            if pdf.startswith('/'):
                pdf = domain + pdf
            pdf_urls.append(pdf)
    return(pdf_urls)

        
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
        link = check_url(data[key])
        pdfs = pdf_anchor_finder(link)
        for pdf in pdfs:
            try:
                pdf_file = urllib.urlopen(pdf).read()
                fp = os.path.join(dir, str(key)) + '.pdf'
                fp_i = 1
                while(os.path.exists(fp)):
                    fp = os.path.join(dir, "{0}_{1}.pdf".format(str(key), fp_i))
                    fp_i = fp_i + 1
                with open(fp, 'wb') as f:
                    f.write(pdf_file)
                pdfs_found += 1
                break
            except Exception as e:
                print(e, key, pdf)
    print('PDFS FOUND: ' + str(pdfs_found))
    

def pdf_finder(data, dir):
    """
    Tries to find as many PDFs as possible from URLs in the dataframe
    :param data: DataFrame
    :param dir: Out directory
    """
    pdf_links = dict(zip(data['i'], data['PDF link']))
    urls = dict(zip(data['i'], data['URL']))
    pdfs_found = 0
    for key in urls:
        fp = os.path.join(dir, str(key)) + '.pdf'
        tell_me = "\033[1m{0}:\033[0m".format(key)
        if os.path.exists(fp):
            tell_me = "{0}\n{1} exists".format(tell_me, fp)
        else:
            link = check_url(pdf_links[key])
            link = link if link else check_url(urls[key])
            if link:
                tell_me = "{0}\noriginal link: {1}".format(tell_me, link)
                pdfs = pdf_anchor_finder(link)
                if len(pdfs):
                    for pdf in pdfs:
                        try:
                            pdf_file = urllib.urlopen(pdf).read()
                            fp_i = 1
                            while(os.path.exists(fp)):
                                fp = os.path.join(dir, "{0}_{1}.pdf".format(str(key), fp_i))
                                fp_i = fp_i + 1
                            with open(fp, 'wb') as f:
                                f.write(pdf_file)
                            try:
                                x = PdfReader(pdf_file)
                                if(len(x.pages) > 1):
                                    tell_me = "{0}\n{1} saved as {2}".format(tell_me, pdf, fp)
                                    pdfs_found += 1
                                else:
                                    tell_me = "{0}\n{1} only one page, not saved".format(tell_me, pdf)
                                    os.remove(pdf_file)
                            except:
                                tell_me = "{0}\n{1} failed to save as valid pdf".format(tell_me, pdf)
                                os.remove(pdf_file)
                        except Exception as e:
                            tell_me = "{0}\n{1}: {2} {3}".format(tell_me, e, pdf, link)
                            continue
            else:
                tell_me = "{0}\n{1}".format(tell_me, "no URL")
        print(tell_me)
        

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
            url = url[26:] if url.startswith(
                'http://scholar.google.com/') else url
            print(i, url)
            unlinkables.append((i, url))
    return unlinkables

