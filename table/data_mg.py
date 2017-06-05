import manager as mg, re
data = mg.get_data()

data = data[data['Sets'] == 'INDI']

paragraphs = mg.get_paragraphs()

for i, row in data.iterrows():
    print(i, row['Title'], row['Authors'], row['Journal'], row['URL'])
    print('\n'.join(paragraphs[i]))
    print('\n\n\n\n\n')

"""
print(scholar_reader.authors(data, 'Sets', split=';'))
data = mg.get_data()
data['Authors'] = data['Authors'].fillna('').apply(lambda x: re.sub('[^ ,&\w]', '', x))
data['Year'] = data['Year'].fillna(0).apply(lambda x: str(int(x))).replace(to_replace='0', value='')
data = data[data['Data Use'] == 'Y']
data = data[data['Journal Category'] == 'Thesis']
table_data = data[['Title', 'Authors', 'Year', 'Citations', 'URL']]
for row in data.iterrows():
    row = row[1]
    print(row['i'], row['Title'], row['URL'], row['Data Use'], row['Sets'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)

invalid_theses = [500, 808, 1038, 1072, 1363]

"""