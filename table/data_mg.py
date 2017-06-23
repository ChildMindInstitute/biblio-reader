import manager as mg
data = mg.get_data()
data = data[data['Data Use'] == 'Y'][data['Journal Category'] == 'Thesis']
table_data = data[['Title', 'Authors', 'Year', 'URL']].sort_values(['Year', 'Title'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)
