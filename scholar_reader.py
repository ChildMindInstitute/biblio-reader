import pandas as pd

with open('FCP_2015.csv', 'r') as f:
    data = pd.read_csv(f)

journal = data['Journal'].dropna().apply(lambda x: str(x).casefold())

year = data['Year'].dropna()

print(journal.value_counts().head(20))

journals = data.groupby('Journal').groups


