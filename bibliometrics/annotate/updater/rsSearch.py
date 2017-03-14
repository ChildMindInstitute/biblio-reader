from Bio import Entrez
import sets, sys
from pybtex.database.input import bibtex
import MySQLdb as mysql

def is_int(a):
	try:
		int (a)
		return True
	except:
		return False

## setup access to pubmed
Entrez.email='drcc@vt.edu'
database='pubmed'

# now we perform the new query
query = """
	(("rest" or "resting" or "resting state" or "spontaneous" or "intrinsic" or
	 "low frequency fluctuation" or "low frequency fluctuations" or
	 "low-frequency fluctuations" or
	 "connectivity" or "connectome" or "default mode" or "\'default mode\'" or
	 "default mode network" or "default network") and ("functional magnetic resonance 
	 imaging" or "functional mri" or "fmri" or "fMRI" or BOLD or "blood oxygenation level dependent"
	 or "blood oxygention level-dependent" or "blood oxygenation-level-dependent" or
	 "blood oxygen level dependent" or "blood oxygen-level dependent" or
	 "blood oxygen-level-dependent" or "functional neuroimaging") or 
	 "resting state functional connectivity MRI" 
	 or 
	 (("resting state" or "rest" or "intrinsic") and ("functional connectivity")) or
	 ((functional or function) and ("connectome" or "connectomes" or "connectomics")) 
	 or 
	 ((("rest" and "state") or "resting state" or "functional connectivity" or 
	 "default mode" or "default network") and ("MRI" or "magnetic resonance imaging")) 
	 or rsfMRI or "resting state fMRI" or "resting state functional magnetic resonance imaging" or
	 "intrinsic connectivity network" or "resting state network" or
	 "resting state networks" or "intrinsic connectivity networks" or
	 "intrinsic brain activity" or "brain functional connectivity")
	 
	 
	 and ("1995/01/01"[Date - Publication] : "2012/12/31"[Date - Publication])
"""

print query

result = Entrez.esearch(db=database,term=(query),retmax=100000)
try:
	record=Entrez.read(result)
except:
	print "Error"

print "\n\nSearch returned %d records"%(int(record['Count']))

new_pmids=set(record['IdList'][:])

too_new=0
not_found=0
rejected=0
already_exist=0
broken = 0
#print "Received %d new pmids"%(len(new_pmids))

# go through and throw out the ones with 
# publication dates > 2012
#
t_new_pmids=new_pmids.copy()
counter = 0
db = mysql.connect(host='localhost', user='root', passwd='renew', db='annotate')
cur = db.cursor()
for id in t_new_pmids:
	cur.execute('select count(*) from rejected_refs where pmid="' + str(id) + '"')
	row = cur.fetchone()
	if int(row[0]) > 0:
		rejected += 1
		new_pmids.remove(id)
		continue	
	cur.execute('select count(*) from refs where pmid="' + str(id) + '"')
	row = cur.fetchone()
	if int(row[0]) > 0:
		already_exist += 1
		new_pmids.remove(id)
		continue
	k_result = Entrez.efetch(db=database, id=id, retmode="xml", \
		rettype="full", tool="bioPython, Bio.Entrez")
	k_record = Entrez.read(k_result, validate=False)

	if k_record:
		try:
			if int(k_record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']) > 2012:
				cur.execute('insert into rejected_refs (pmid) VALUES ("' + str(id) + '")')
				new_pmids.remove(id)
				too_new=too_new+1
		except:
			print "%d/%d Could not find pub date for %s"%(counter, len(t_new_pmids), id)
			cur.execute('insert into rejected_refs (pmid) VALUES ("' + str(id) + '")')
			new_pmids.remove(id)
			broken += 1
	else:
		print "%d/%d PMID %s not found"(counter, len(t_new_pmids), id)
		cur.execute('insert into rejected_refs (pmid) VALUES ("' + str(id) + '")')
		new_pmids.remove(id)
		not_found=not_found+1
	if counter % 100 == 0:
		print "%d/%d Finished"%(counter, len(t_new_pmids))
	counter = counter+1

print "%d were previously rejected\n%d already exist\n%d pmids were not found\n%d had no publication date\n%d were published after 2012\nthere are %d remaining"%(
	rejected, already_exist, not_found, broken, too_new, len(new_pmids))
db.commit()
db.close()

outf = open(sys.argv[1], 'w')
for pmid in new_pmids:
	outf.write(str(pmid) + "\n")

