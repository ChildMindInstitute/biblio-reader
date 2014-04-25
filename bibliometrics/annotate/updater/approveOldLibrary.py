import sys
from collections import *
import pybtex.database.input.bibtex
import MySQLdb as mysql

tagTr = {"Meta-Analysis/Reviews": "review meta analysis",
		"Animal Model": "animal models"}

f = sys.argv[1]
parser = pybtex.database.input.bibtex.Parser()
bib_data = parser.parse_file(f)

db = mysql.connect(host='localhost', user='root', passwd='renew', db='annotate')
cur = db.cursor()

def approve(bib_data, cur):
	for k in bib_data.entries.keys():
		if not 'pmid' in bib_data.entries[k].fields:
			continue
		pmid = bib_data.entries[k].fields['pmid']
		cur.execute('select ref_id from refs where pmid="' + str(pmid) + '" limit 1')
		row = cur.fetchone()
		if row is None or len(row) != 1:
			print "Skipped PMID", pmid
			continue
		refsId = int(row[0])
		query = 'update expert_assignments set status="pending_tagging" where ref_id=' + str(refsId)
		print query
		cur.execute(query)
def tag():
	for k in bib_data.entries.keys():
		if not 'pmid' in bib_data.entries[k].fields:
			continue		
		pmid = bib_data.entries[k].fields['pmid']
		
		if not 'mendeley-tags' in bib_data.entries[k].fields:
			continue
		oldTags = bib_data.entries[k].fields['mendeley-tags']
		newTagList = oldTags.split(",")
		trTagList = []
		for tag in newTagList:
			if tag in tagTr:
				tag = tagTr[tag]
			tag = tag.lower().replace(" ", "_")
			tag = '"' + tag + '"'
			trTagList.append(tag)		
		newTags = ",".join(trTagList)
		
		cur.execute('select ref_id from refs where pmid="' + str(pmid) + '" limit 1')
		row = cur.fetchone()
		if row is None or len(row) != 1:
			print "Skipped PMID", pmid
			continue
		refsId = int(row[0])
		#query = 'update expert_assignments set tags = concat_ws(",", tags, ' + newTags + ') where ref_id=' + str(refsId)
		#print query
		#cur.execute(query)
		#query = 'update expert_assignments set status="tagging_complete" where ref_id=' + str(refsId)
		#print query
		query = 'update refs set mendeley_tags="' + oldTags + '", source="old_library" where ref_id=' + str(refsId)
		cur.execute(query)		
tag()		
db.commit()
db.close()	

