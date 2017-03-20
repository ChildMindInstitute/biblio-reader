import sys
from collections import *
import pybtex.database.input.bibtex

f = sys.argv[1]
parser = pybtex.database.input.bibtex.Parser()
bib_data = parser.parse_file(f)

fields = ['pmid', 'title', 'author', 'journal', 'year', 'volume', 'pages', 'number', 'abstract', 'affiliation', 
		'doi', 'isbn', 'journal_iso', 'mesh_heading_list', 'master_tags']
taggedPhrases = {'support vector machine': 'supervised_learning',
				'svm': 'supervised_learning',
				' monkey ': 'animal_models',
				' rat ': 'animal_models'}
allTags = ['clinical','basic_neuroscience','brain_and_behavior','functional_anatomy','multimodal','genetics','animal_models','methodology','review_meta_analysis','1000_functional_connectomes','graph_theory','supervised_learning','unsupervised_learning','independent_component_analysis','seed_based_correlation','mark_for_removal','expert_needed']
ignoredFields = ['__markedentry']
refs = []

for k in bib_data.entries.keys():
	currentRef = defaultdict(str)
	for field in bib_data.entries[k].persons.keys():
		currentRef[field] = " and ".join([unicode(p) for p in bib_data.entries[k].persons[field]])
	for field in bib_data.entries[k].fields.keys():
		assert field in fields or field in ignoredFields, parts
		currentRef[field] = bib_data.entries[k].fields[field]
	refs.append(currentRef)

i = 0
outf = open(f.split(".")[0] + ".sql", 'w')
for ref in refs:
	i += 1
	tags = set()
	for tag in allTags:
		if ' ' + tag.replace('_', ' ') + ' ' in ref['abstract']:
			tags.add(tag)
	ref['master_tags'] = ','.join(tags) 
	
	line = []
	line.append('"' + str(i) + '"')
	for field in fields:
		line.append('"' + ref[field].strip() + '"')
	outf.write("INSERT INTO refs")
	outf.write(" (" + ", ".join([('`' + field + '`') for field in fields]) + ")")
	outf.write(" VALUES (" + ", ".join([('"' + ref[field].strip().encode('utf8') + '"') for field in fields]) + ");\n")
	#print ",".join(line)	
