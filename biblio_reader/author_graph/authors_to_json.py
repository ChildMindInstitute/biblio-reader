import manager as mg, json, os
data = mg.get_data()
from biblio_reader import scholar_reader

types = {auth: 'Contributor' in val for auth, val in scholar_reader.authors(data, 'Contributor').items()}

with open(os.path.join(mg.dir(os.path.join('data', 'author-links')), 'objects.json'), 'w') as o:
    json.dump([{'name': auth,
                'type': types[auth],
                'depends': list({aff for aff in affils if aff != auth and aff != 'others'})}
               for auth, affils in scholar_reader.authors(data, 'Authors', split=' & ').items() if auth != 'others'],
              o, sort_keys=True, indent=4)




