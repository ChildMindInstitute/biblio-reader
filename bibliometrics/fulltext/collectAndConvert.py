#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collectAndConvert.py

Script to collect the pdfs that we can from the previously conducted Google
    Scholar search and convert them all to text files that we can work with.

Requires convertToText.py by Matt Doherty which in turn requires pdftotext.

Author:
	– Jon Clucas, 2017 (jon.clucas@childmind.org)
"""

import convertToText, os, pandas as pd, re, urllib

def main():
    # make in and out dirs if they don't already exist
    for subdir in ['pdf', 'txt']:
        if not os.path.exists(os.path.join('fulltexts', subdir)):
            os.makedirs(os.path.join('fulltexts', subdir))
    
    # get the list of pdf files to work with
    with open('../../FCP_2015.csv', 'r', encoding='utf-8') as f:
        data = pd.read_csv(f)

    # get the pdfs
    for pdf, authors, title in zip(data['PDF link'], data['Authors'], data[
                               'Title']):
        if isinstance(pdf, str):
            pdf = pdf.replace('http://scholar.google.com/', '')
            file = os.path.join('fulltexts', 'pdf', ''.join([re.sub(r'\W', '',
                   authors.split(' & ')[0])[:8], '_', re.sub(r'\W', '', title)[
                   :8], '.pdf']))
            if not os.path.exists(file) or os.path.getsize(file) == 0:
                try:
                    pdffile = urllib.request.urlopen(pdf)
                    print(file)
                    with open(file, 'wb') as ofile:
                        ofile.write(pdffile.read())
                except Exception as e:
                    print('\t'.join([str(e), pdf]))
                    next
                
    
    # convert the pdfs to txts
    convertToText.walkAndText(os.path.abspath(os.path.join('fulltexts', 'pdf')
                             ), os.path.abspath(os.path.join('fulltexts',
                             'txt')))
		
# ============================================================================
if __name__ == '__main__':
    main()