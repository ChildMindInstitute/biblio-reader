import pandas as pd
import matplotlib.pyplot as plt
import manager, os, datetime, re
data = manager.get_data()


def countstats(series, path=None, row_limit=None):
    stat = series.dropna().apply(lambda x: str(x).casefold()).value_counts()
    if row_limit:
        stat = stat.head(row_limit)
    if path:
        stat.to_csv(path=path, header=True)
    else:
        return stat


def citations_per_year(data, sort=False):
    data['CPY'] = data['Citations'] / (datetime.datetime.now().year + 1 - data['Year'])
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


def count_item_year(data, column, item):
    return data.loc[data[column].apply(lambda x:str(x).casefold()) == item]['Year'].value_counts()


def year_sum_graph(data, column):
    df = data.groupby('Year').sum()[column]
    plt.figure()
    plt.ylabel(column)
    plt.xlabel('Year')
    df.plot.bar()
    plt.savefig('Citations_sum.png', bbox_inches='tight')


def categorize_journals(data, categories):
    res = {}
    if len(os.listdir(categories)) == 0:
        print('No journal categories')
        return
    for category in os.listdir(categories):
        cat_name = category.replace('.txt', '')
        with open(os.path.join(categories, category)) as c:
            keywords = [keyword.strip() for keyword in c.readlines()]
        if '_URL' in cat_name:
            type = 'URL'
            cat_name = cat_name.replace('_URL', '')
        else:
            type = 'Journal'
        data_journals = data.dropna(subset=[type])
        for keyword in keywords:
            res.update({row[1]['i']: cat_name for row in data_journals.iterrows() if keyword in row[1][type].lower()})
    res.update({i: 'Unknown' for i in range(0, len(data)) if i not in res})
    return res
