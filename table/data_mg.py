import manager as mg, re

data = mg.get_data()

"""
data = data[data['Data Use'] == 'Y'][data['Journal Category'] == 'Thesis']
table_data = data[['Title', 'Authors', 'Year', 'URL']].sort_values(['Year', 'Title'])
table_data.to_csv(path_or_buf='Table_data.csv', index=False)
"""

x = re.sub('[ \t]+', ',', """1	Y
4	Y
5	Y
19	Y
23	Y
25	Y
27	Y
29	Y
30	Y
35	Y
39	Y
43	Y
50	Y
57	Y
68	Y
78	Y
81	Y
87	Y
97	Y
116	Y
124	Y
135	Y
149	Y
179	N
183	Y
192	Y
193	Y
265	Y
267	Y
275	Y
320	Y
330	Y
351	Y
392	Y
425	Y
460	Y
514	Y
522	Y
629	Y
645	Y
695	Y
857	Y
1009	Y
1020	Y
1063	Y
1093	I
1128	Y
1225	Y
1240	Y
1333	Y
1348	I
1353	I
1366	I
1438	Y
1455	Y
1465	Y
1490	I""").strip()

print(data.loc[1, 'Data Use'])

"""
x = [tuple(r.split(',')) for r in x.split('\n')]
for i, use in x:
    data.loc[i, 'Data Use'] = use
"""
data['i'] = data['i'].apply(lambda x: int(x))
mg.update_data()
