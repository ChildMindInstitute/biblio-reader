from Bio import Entrez
from collections import defaultdict
import re, sys
import MySQLdb as mysql
#from sets import set
	 
db = mysql.connect(host='127.0.0.1', user='root', passwd='renew', db='annotate')
cur = db.cursor()

Entrez.email='drcc@vt.edu'
database='pubmed'
fmriQuery = """
("fMRI" or "functional magnetic resonance imaging") and ("YEAR_HERE"[Publication Date])
"""

print "\t".join(["year", "fmriTotal", "rsTotal"])
for year in range(1995, 2013):
	cur.execute("select count(*) from (select refs.pmid from expert_assignments, refs where status='tagging_complete' and year='%s' and refs.ref_id=expert_assignments.ref_id group by refs.pmid) as x" % (year))
	row = cur.fetchone()
	if row is None:
		rsTotal = 0
	else:
		rsTotal = int(row[0])	
	
	query = fmriQuery.replace("YEAR_HERE", str(year))
	searchResult = Entrez.esearch(db = database, term = query)
	searchRecord = Entrez.read(searchResult)

	# get the number of matches to the query
	fmriTotal = searchRecord['Count']

	print "\t".join([str(x) for x in year, fmriTotal, rsTotal])

