import sets
import pybtex.database.input.bibtex
import pybtex.database.output.bibtex

## setup bibtex parser
writer = pybtex.database.output.bibtex.Writer()
# write_file(self, bib_data, filename)

print "Reading done"
parser = pybtex.database.input.bibtex.Parser()
bib_done=parser.parse_file('done_20130709.bib')

print "Reading original"
parser = pybtex.database.input.bibtex.Parser()
bib_orig=parser.parse_file('new_pmids.bib')

lib_pmids=set()

new_entries = {}
for k in bib_done.entries.keys():
	if k in bib_orig.entries:
		new_entries[k] = bib_orig.entries[k]
		
bib_done.entries = new_entries

print "Writing rebib done"
writer.write_file(bib_done, 'rebib_done_20130709.bib')