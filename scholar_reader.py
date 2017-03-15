import pandas as pd

with open('FCP_2015.csv', 'r') as f:
    data = pd.read_csv(f)


def countstats(series, path=None, row_limit=None):
    stat = series.dropna().apply(lambda x: str(x).casefold()).value_counts()
    if row_limit:
        stat = stat.head(row_limit)
    if path:
        stat.to_csv(path=path, header=True)
    else:
        print(stat)


def citations_per_year(sort=False):
    data['CPY'] = data['Citations'] / (2018 - data['Year'])
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)


def count_tuples():
    journal_and_year = pd.Series(list(zip(data['Journal'].dropna(), data['Year'].dropna())))
    countstats(journal_and_year)
count_tuples()



