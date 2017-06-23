# biblio_reader

biblio_reader is the main analysis utility and provides support from the following modules:
## scholar_reader.py
Utilities to analyze trends in the main CSV file including the following:
* Matplotlib graph tools for counting pandas columns and visualizing data
* Updates to the CSV file with average no. of citations per year for each publication
* Journal attribute counter for each journal article
* Author relationship dictionary generator
* Statistics and impacts printer

These utilities are independent and ?  useful throughout the project

## pdf_utils.py

Utilities to help find full-text PDFs from each publication URL
 
 Unfortunately this cannot search every publication programatically as some URLs are protected by web-crawling blocks  

## text_tools

* Converts PDFs to texts
* Finds specific paragraphs in each full text that contain the keywords Google Scholar matched
* Associates each article with the keyword sets the user inputs in manager.py

## article_assign.py
Assigns articles randomly to be reviewed manually by multiple reviewers
## review_analysis.py
Analyzes the reviewer inputs for each article, and categorizes publications accordingly
## pub_med.py
Finds affiliation and topic information for articles that are able to be found on Pubmed Central
## map
Takes affiliation information parsed on pub_med.py and maps it using [Google Maps Javascript API](https://developers.google.com/maps/documentation/javascript/)