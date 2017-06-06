import manager as mg, re
data = mg.get_data()


data = mg.get_data()
data = data[data['Data Use'] == 'Y']
data = data[data['Journal Category'] == 'Thesis']
table_data = data[['Title', 'Authors', 'Year', 'Citations', 'URL']]
for row in data.iterrows():
    row = row[1]
    print(row['i'], row['Title'], row['URL'], row['Data Use'], row['Sets'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)

invalid_theses = [500, 808, 1038, 1072, 1363]

"""