import PyPDF2, os, re
import pandas as pd
import urllib.parse as urlparse
import urllib.request as urllib
import manager as mg
from bs4 import BeautifulSoup as bs
from biblio_reader.text_tools import convertToText
data = mg.get_data()

def pdfopener(data, dir):
    for row in data.iterrows():
        url = row[1]['URL']
        if pd.isnull(url):
            continue
        if '.pdf' in url:
            try:
                pdf = urllib.urlopen(url)
                with open(dir + str(row[1]['i']) + '.pdf', 'wb') as f:
                    f.write(pdf.read())
            except:
                print('Unable to open:', row[1]['i'])
                print(url)
                continue


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


def find_corrupted(pdf_directory):
    res = []
    for path, dirs, files in os.walk(pdf_directory):
        for file in files:
            full_file = '/'.join([path, file])
            try:
                PyPDF2.PdfFileReader(full_file)
            except:
                res.append(int(file.replace('.pdf', '')))
    return res


def find_paragraphs(txt_directory, terms, outfile=None):
    res = {}
    terms = list(map((lambda x: x.lower()), terms))
    for path, dirs, files in os.walk(txt_directory):
        for file in files:
            full_file = '/'.join([path, file])
            with open(full_file, 'r') as f:
                text = f.read()
            paragraphs = re.split(r'[ \t\r\f\v]*\n[ \t\r\f\v]*\n[ \t\r\f\v]*', text)
            paragraphs = [paragraph.lower() for paragraph in paragraphs if isinstance(paragraph, str)]
            key_paragraphs = []
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs += [paragraph for paragraph in paragraphs
                                   if re.search(re_term, paragraph) and paragraph not in key_paragraphs]
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs = [str(re.sub(re_term, '@@@@' + term, paragraph)) for paragraph in key_paragraphs]
            res[int(file.replace('.txt', ''))] = key_paragraphs
    if outfile:
        with open(outfile, 'w') as f:
            f.write(str(res))
    return res

def main():

    PDF_DIR = mg.dir(os.path.join(mg.INPUT_PATH, 'pdfs'))

    TXT_DIR = mg.dir(os.path.join(mg.ROOT_PATH, 'txts'))

    if len(os.listdir(PDF_DIR)) == 0:
        pass

    dict_data = dict(zip(data['i'], data['URL']))
    dict_titles = dict(zip(data['i'], zip(data['Title'], data['URL'])))
    valid_data = {key: value for key, value in dict_data.items() if not isinstance(value, float)}
    pdfs = [int(pdf.replace('.pdf', '')) for pdf in os.listdir('../inputs/pdfs') if not pdf.startswith('.')]
    no_pdfs = {key: value for key, value in dict_data.items() if key not in pdfs}
    print(len(pdfs))
    convertToText.walkAndText(mg.dir(os.path.join(mg.INPUT_PATH, 'pdfs')), '../outputs/txts')
    paragraph_dict = find_paragraphs('../outputs/txts', mg.FCP_TERMS, outfile='../outputs/paragraphs.txt')
    empty_paragraphs = [key for key, value in paragraph_dict.items() if len(value) == 0]
    bad_data = list(no_pdfs) + empty_paragraphs
    data[data['i'].isin(bad_data)].to_csv(path_or_buf='../outputs/unlinkables.csv', index=False)

