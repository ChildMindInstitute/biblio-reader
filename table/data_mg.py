import manager as mg
import collections

data = mg.get_data()

data = data[data['Data Use'] != 'I']
table_data = data[['Title', 'Authors', 'Journal', 'Year', 'CPY', 'URL']]

table_data.to_csv(path_or_buf='Table_data.csv', index=False)
