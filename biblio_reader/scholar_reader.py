import pandas as pd
import matplotlib.pyplot as plt
import manager as mg
import os, datetime, collections, ast
data = mg.get_data()


def count_visualizer(value_count, out, stat_type='csv', row_limit=None):
    """
    Counts values of specific columns in dataframe
    :param value_count: A value count series (see pandas value_count function)
    :param out: output file name
    :param stat_type: one of: csv, bar, pie
    :param row_limit: Sets a limit to how many highest values should be counted
    :return: csv, bar, or pie file
    """
    if row_limit:
        value_count = value_count.head(row_limit)
    if stat_type == 'csv':
        value_count.to_csv(path=out, header=[value_count.name, 'Count'])
        return
    plt.figure()
    plt.ylabel('Count')
    plt.title(value_count.name + ' Count')
    if stat_type in ['pie',  'bar']:
        value_count.plot(kind=stat_type)
        plt.xlabel('', visible=False)
    else:
        raise IOError("Invalid stat type")
    plt.savefig(out, bbox_inches='tight')
    plt.clf()




def citations_per_year(data, sort=False):
    data['CPY'] = data['Citations'] / (datetime.datetime.now().year + 1 - data['Year'])
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)


def value_count_graph(data, column, xcount=10):
    series = data[column]
    value_series = value_counter(data, series)
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


def count_sets(data):
    if 'Sets' not in data:
        raise IOError('No sets to work with')
    sets = [dset for dsets in data['Sets'].dropna() for dset in dsets.split(';')]
    return pd.Series(collections.Counter(sets), name='Sets')


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

#countstats(data['Journal'], path=os.path.join(manager.OUTPUT_PATH, 'stats/journalstat.csv'))

#count_visualizer(count_sets(data), , row_limit=20, stat_type='pie')
"""

journal_dict = {journal.lower(): typ for journal, typ in journal_dict.items()}
with open(os.path.join(manager.OUTPUT_PATH, 'stats/journal_use_stat.csv'), 'r') as f:
    journal_stat = pd.read_csv(f)
    res = []
    for journal in journal_stat['Journal']:
        try:
            res.append(journal_dict[journal])
        except KeyError:
            res.append('Unknown')
    journal_stat['Journal Type'] = res
    journal_stat.to_csv(path_or_buf=os.path.join(manager.OUTPUT_PATH, 'stats/journal_use_stat.csv'), index=False)
"""


