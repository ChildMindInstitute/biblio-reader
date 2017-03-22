I got tf–idf working from [Cameron](https://github.com/mfleisch5/scholar.py/commits/master/bibliometrics?author=ccraddock) & [Matt](https://github.com/mfleisch5/scholar.py/commits/master/bibliometrics?author=mattdoherty)'s [code](https://github.com/ccraddock/restingstate_bibliometrics), although there's plenty left for you to do with it (like tweak the terms, outputs, and visualizations).

To get it to work takes a little pre-work:

- poppler for OS X includes pdftotext, which you need to get the text from the pdfs
- I had to run 2to3 on textmining before it would work since it was built in Python2 and hasn't had a release yet for Python3

Once you do all that, you have to run `scholar.py/bibliometrics/fulltext/collectAndConvert.py` to get the pdfs and their full texts before you run `scholar.py/bibliometrics/tfidf/tfidf.py`.
