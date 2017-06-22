import manager as mg, json, os
data = mg.get_data()
data = data[data['Data Use'] == 'Y'].dropna(subset=['Authors'])
from biblio_reader import scholar_reader


types = {auth: 'Contributor' in val for auth, val in scholar_reader.authors(data, 'Contributor').items()}
contributors = {author for i, authors in zip(data['i'], data['Authors']) for author in authors.split(' & ')
                            if i in {1, 5, 74, 92, 653}}
with open(os.path.join(mg.dir(os.path.join('data', 'author-links')), 'objects.json'), 'w') as o:
    json.dump([{'name': auth,
                'type': "Contributor" if auth in contributors else 'Not a Contributor',
                'depends': list({aff for aff in affils if aff != auth and aff != 'others'})}
               for auth, affils in scholar_reader.authors(data, 'Authors', split=' & ').items() if auth != 'others'],
              o, sort_keys=True, indent=4)
