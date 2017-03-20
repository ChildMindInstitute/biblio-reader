import pandas as pd
import matplotlib.pyplot as plt


with open('FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)


def countstats(series, path=None, row_limit=None):
    stat = series.dropna().apply(lambda x: str(x).casefold()).value_counts()
    if row_limit:
        stat = stat.head(row_limit)
    if path:
        stat.to_csv(path=path, header=True)
    else:
        return stat


def citations_per_year(sort=False):
    data['CPY'] = data['Citations'] / (2018 - data['Year'])
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)

def value_count_graph(series, xcount=10):
    value_series = countstats(series)
    value_df = pd.DataFrame()
    for item in value_series.head(xcount).index.tolist():
        value_df[item] = count_item_year(series.name, item)
    value_df.sort_index(inplace=True, ascending=False)
    value_df = value_df.T
    plt.figure()
    plt.xlabel(series.name)
    plt.ylabel('Number of' + series.name + 's')
    value_df.plot.bar(stacked=True)
    plt.savefig('Journals_by_year.png', bbox_inches='tight')


def count_item_year(column, item):
    return data.loc[data[column].apply(lambda x:str(x).casefold()) == item]['Year'].value_counts()


def year_sum_graph(column):
    df = data.groupby('Year').sum()[column]
    plt.figure()
    plt.ylabel(column)
    plt.xlabel('Year')
    df.plot.bar()
    plt.savefig('Citations_sum.png', bbox_inches='tight')

value_count_graph(data['Journal'], xcount=20)
