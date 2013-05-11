from Bio import Entrez
import sets
from pybtex.database.input import bibtex

def is_int(a):
    try:
        int (a)
        return True
    except:
        return False

## setup access to pubmed
Entrez.email='drcc@vt.edu'
database='pubmed'

## setup bibtex parser
parser = bibtex.Parser()

bib_data=parser.parse_file('cmi_library.bib')

id_num=0
lib_pmids=set()

for k in bib_data.entries.keys():
    if 'pmid' in bib_data.entries[k].fields.keys():
        lib_pmids.add(bib_data.entries[k].fields['pmid'])
    else:
        print "PubMed ID missing for %s: %s"%(k, \
            bib_data.entries[k].fields['title'])

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

new_pmids=set(record['IdList'])

too_new=0
not_found=0
print "Received %d new pmids"%(len(new_pmids))

# go through and throw out the ones with 
# publication dates > 2012
#
t_new_pmids=new_pmids.copy()
counter = 0
for id in t_new_pmids:
    k_result = Entrez.efetch(db=database, id=id, retmode="xml", \
        rettype="full", tool="bioPython, Bio.Entrez")
    k_record = Entrez.read(k_result)

    if k_record:
        try:
            if int(k_record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']) > 2012:
                new_pmids.remove(id)
                too_new=too_new+1
        except:
            print "%d/%d Could not find pub date for %s"%(counter, len(t_new_pmids), id)
    else:
        print "%d/%d PMD %s not found"(counter, len(t_new_pmids), id)
        new_pmids.remove(id)
        not_found=not_found+1
    if counter % 100 == 0:
        print "%d/%d Finished"%(counter, len(t_new_pmids))
    counter = counter+1

print "%d pmids were not found, %d were published after 2012, there are %d remaining"%(
    not_found, too_new, len(new_pmids))

with open('new_pmids.txt', 'w') as f:
    f.write("\n".join(new_pmids))

new_diff_pmids=new_pmids.difference(lib_pmids)

print "%d of which are not in the old library"%(len(new_diff_pmids))
with open('new_pmids_diff.txt', 'w') as f:
    f.write("\n".join(new_diff_pmids))

old_diff_pmids=lib_pmids.difference(new_pmids)
print "checking %d pmids in the old library that are not in the new list"% \
    (len(old_diff_pmids))

t_old_diff_pmids=old_diff_pmids.copy()
for id in t_old_diff_pmids:
    k_result = Entrez.efetch(db=database, id=id, retmode="xml", \
        rettype="full", tool="bioPython, Bio.Entrez")
    k_record = Entrez.read(k_result)

    if k_record:
        try:
            if int(k_record[0]['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']) > 2012:
                old_diff_pmids.remove(id)
                too_new=too_new+1
        except:
            print "Could not find pub date for %s"%(id)
    else:
        old_diff_pmids.remove(id)
        not_found=not_found+1

print "of %d pmids in the old library not in the new list, %d were too new, and %d have pmids that could not be found"%(len(old_diff_pmids), too_new, not_found)
with open('old_pmids_diff.txt', 'w') as f:
    f.write("\n".join(old_diff_pmids))


