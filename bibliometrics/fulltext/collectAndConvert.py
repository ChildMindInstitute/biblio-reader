import convertToText, os, pandas as pd, re, urllib

def main():
    with open('../../FCP_2015.csv', 'r', encoding='utf-8') as f:
        data = pd.read_csv(f)
		
    for pdf, authors, title in zip(data['PDF link'], data['Authors'], data[
                               'Title']):
        if isinstance(pdf, str):
            pdf = pdf.replace('http://scholar.google.com/', '')
            file = os.path.join('fulltexts', ''.join([re.sub(r'\W', '',
                   authors.split(' & ')[0])[:8], '_', re.sub(r'\W', '', title)[
                   :8], '.pdf']))
            if not os.path.exists(file):
                try:
                    pdffile = urllib.request.urlopen(pdf)
                except Exception as e:
                    print(e)
            
                print(file)
                with open(file, 'wb') as ofile:
                    ofile.write(pdffile.read())
		
# ============================================================================
if __name__ == '__main__':
    main()