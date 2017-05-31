import pandas as pd
import matplotlib.pyplot as plt
import manager as mg
import os, datetime, collections, re
data = mg.get_data()


def count_visualizer(value_count, out, stat_type, name, row_limit=None):
    """
    Counts values of specific columns in dataframe
    :param value_count: A value count series, dict, or LOT (see pandas value_count function)
    :param out: output file name
    :param stat_type: one of: csv, bar, pie
    :param row_limit: Sets a limit to how many highest values should be counted
    :return: csv, bar, or pie file
    """
    value_count = dict(value_count)
    if row_limit:
        value_count = sorted(value_count, key=value_count.get)[:row_limit]
    plt.figure()
    if stat_type == 'bar':
        plt.bar(range(len(value_count)), list(value_count.values()), align='center')
        plt.xticks(range(len(value_count)), value_count.keys())
    elif stat_type == 'pie':
        plt.pie(list(value_count.values()), labels=value_count.keys(), autopct='%1.1f%%', shadow=True)
        plt.axis('equal')
    else:
        raise IOError('Invalid stat type')
    plt.title(name)
    plt.savefig(out + '.png', bbox_inches='tight')
    plt.clf()


def citations_per_year(data, sort=False):
    data['CPY'] = data['Citations'] / (datetime.datetime.now().year + 1 - data['Year'])
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)

def stack_bar(data, column, stacker, out, split=None):
    plt.figure()
    last_stack = None
    for stack in data[stacker].value_counts().index.tolist():
        print(stack)
        if not split:
            stacked = [value for value in data[data[stacker] == stack].dropna()[column]]
        else:
            stacked = [value for values in
                       data[data[stacker] == stack].dropna()[column]
                       for value in values.split(split)]
        print()
        stacked = collections.Counter(stacked)
        print(*stacked.items(), sep='\n')
        plt.bar(range(len(stacked)), list(stacked.values()), align='center', label=stack, bottom=last_stack)
        #plt.set_prop_cycle(cycler('color', ['r', 'g', 'b', 'y']))
        last_stack = stack
    plt.xticks(range(len(stacked)), stacked.keys())
    plt.savefig(out)

#stack_bar(data, 'Sets', 'Year', 'tester.png', split=';')


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

def author_links(data):
    links = {}
    data = data.dropna(subset=['Authors'])
    for i, authors in zip(data['i'], data['Authors']):
        authors = [re.sub('[^\w ,]', '', author) for author in authors.split(' & ') if author != 'others']
        for author in authors:
            others = {auth for auth in authors if auth != author}
            if author in links:
                links[author].update(others)
            else:
                links[author] = others
    return links


author_counts = {}

#countstats(data['Journal'], path=os.path.join(manager.OUTPUT_PATH, 'stats/journalstat.csv'))

#count_visualizer(count_sets(data), 'tester', 'pie', 'Set Use')

#print(list(count_sets(data).index))
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