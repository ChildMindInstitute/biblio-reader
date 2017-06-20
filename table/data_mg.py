import manager as mg, json, re

data = mg.get_data()

"""
data = data[data['Data Use'] == 'Y'][data['Journal Category'] == 'Thesis']
table_data = data[['Title', 'Authors', 'Year', 'URL']].sort_values(['Year', 'Title'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)
"""

with open(mg.ROOT_PATH + '/x.txt', 'r') as f, open(mg.ROOT_PATH + '/z.txt', 'r') as j, \
    open(mg.ROOT_PATH + '/y.txt', 'r') as jourls:
    x = [r.strip() for r in f.readlines()]
    r = [int(l.strip()) for l in j.readlines()]
    journals = [j.strip() for j in jourls.readlines()]
y = dict(zip(r, x))



with open(mg.ROOT_PATH + '/journals.csv', 'r') as lf:
    l = [tuple([journ] + xl.rstrip('; ,,\n').split(',')) for journ, xl in zip(journals, lf.readlines())]


l = [x for x in l if len(x) == 3]
l = {journal: {'CiteScore': impact,
               'Categories': [y[int(cat)] for cat in categories.split(';')]}
     for journal, impact, categories in l}




with open(mg.INPUT_PATH + '/journal_impacts.json', 'w') as js:
    json.dump(l, js, indent=4)