import pandas as pd
import matplotlib.pyplot as plt
import os, csv


with open('../inputs/FCP_DATA.csv', 'r') as f:
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


def checker_directory(directory):
    checks = {}
    for check in os.listdir(directory):
        full_path = '/'.join([directory, check])
        with open(full_path, 'r') as f:
            reader = list(csv.reader(f))
            for rows in reader[1:]:
                k = int(rows[0])
                v = (check, rows[1].replace(' and ', '').upper())
                if k not in checks:
                    checks[k] = [v]
                else:
                    checks[k].append(v)
    return checks

"""
checker = checker_directory('../inputs/Article_Checks')
double_checked = {key: check for key, check in checker.items() if len(check) > 1}
"""
print(len([year for year in data['Year'] if year == 2016]))
