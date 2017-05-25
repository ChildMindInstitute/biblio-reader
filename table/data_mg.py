import manager as mg
import collections

data = mg.get_data()

table_data = data[['Title', 'Authors', 'Journal', 'Year', 'CPY', 'URL']].copy()

table_data.to_csv(path_or_buf='Table_data.csv', index=False)
