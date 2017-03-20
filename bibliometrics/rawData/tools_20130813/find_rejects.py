import sets
from pybtex.database.input import bibtex

## setup bibtex parser
parser = bibtex.Parser()
bib_current = parser.parse_file('rebib_done_20130709.bib')

parser = bibtex.Parser()
bib_new = parser.parse_file('new_pmids.bib')

current_pmids = set([bib_current.entries[k].fields['pmid'] for k in bib_current.entries.keys()])

outf = open('output.sql', 'w')
for k in bib_new.entries.keys():
	if bib_new.entries[k].fields['pmid'] not in current_pmids:
		if int(bib_new.entries[k].fields['year']) <= 2007:
			outf.write("INSERT INTO rejected_refs (`pmid`) VALUES (" + str(bib_new.entries[k].fields['pmid']) + ");\n")
        