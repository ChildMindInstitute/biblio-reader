import os
import sys
import json
br_path = os.path.abspath(os.path.join(__file__,
                          os.pardir, os.pardir))
if br_path not in sys.path:
    sys.path.append(br_path)
import manager as mg
data = mg.get_data()
data = data[(data['Data Use'] == 'Y')].dropna(subset=['Authors'])
from biblio_reader import scholar_reader


types = {auth: 'Contributor' in val for auth, val in scholar_reader.authors(data, 'Contributor').items()}
contributors = {author for i, authors in zip(data['i'], data['Authors']) for author in authors.split(' & ')
                            if i in mg.CONTR_PAPERS}


with open(os.path.join(mg.dir(os.path.join('data', 'author-links')), 'objects.json'), 'w') as o:
    json.dump([{'name': auth,
                'type': "Contributor" if types[auth] is True else 'Not a Contributor',
                'depends': list({aff for aff in affils if aff != auth and aff != 'others'})}
               for auth, affils in scholar_reader.authors(data, 'Authors', split=' & ').items() if auth != 'others'],
              o, sort_keys=True, indent=4)
