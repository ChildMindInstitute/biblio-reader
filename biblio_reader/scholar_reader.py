import pandas as pd
import matplotlib.pyplot as plt
import manager as mg
import os, datetime, collections, re
import numpy as np
data = mg.get_data()
data = data[data['Data Use'] != 'I']


def count_visualizer(value_count, stat_type, name, out, row_limit=None):
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


def stack_bar(data, column, stacker, stack_type, out, split=None):
    plt.figure()
    stacks = []
    for stack in stacker:
        if not split:
            stacked = [value for value in data[data[stack_type] == stack][column].dropna()]
        else:
            stacked = [value for values in
                       data[data[stack_type] == stack][column].dropna()
                       for value in values.split(split)]
        stacks.append({typ: count for typ, count in sorted(collections.Counter(stacked).items())[:10]})
    max_stack = max(stacks, key=len)
    last_stack = list(np.zeros(len(max_stack), dtype=np.int))
    for name, stack in zip(stacker, stacks):
        empty_stacks = {key: 0 for key in max_stack.keys() if key not in stack}
        if len(empty_stacks) > 0:
            stack.update(empty_stacks)
            stack = {key: value for key, value in sorted(stack.items())}
        print(name, stack, last_stack)
        plt.bar(range(len(stack)), list(stack.values()), align='center', label=name,
                bottom=last_stack)
        last_stack = [x + y for x, y in zip(last_stack, list(stack.values()))]
    plt.xticks(range(len(max_stack)), max_stack.keys())
    plt.title(column + ' by ' + stack_type)
    plt.legend()
    plt.savefig(out)

    """
        if len(last_stack) > len(stacked):
            stacked.update({'No data ' + str(i): 0 for i in range(len(last_stack) - len(stacked))})
        last_stack += [0 for i in range(len(stacked) - len(last_stack))]
        print(len(stacked), len(last_stack))
        plt.bar(range(len(stacked)), list(stacked.values()), align='center', label=stack,
                bottom=last_stack)
        last_stack = [x + y for x, y in zip(last_stack, list(stacked.values()))]
    plt.xticks(range(len(stacked)), stacked.keys())
    plt.title(column + ' by ' + stack_type)
    plt.legend()
    plt.savefig(out)
    """

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

STAT_DIR = mg.dir(os.path.join(mg.OUTPUT_PATH, 'stats'))
stack_bar(data, 'Data Use', range(2009, 2018), 'Year', os.path.join(STAT_DIR, 'sets_by_year.png'))
"""
STAT_DIR = mg.dir(os.path.join(mg.OUTPUT_PATH, 'stats'))
count_visualizer(data['Year'].dropna().apply(lambda x: int(x)).value_counts(), 'bar', 'Year Count',
                 os.path.join(STAT_DIR, 'Year_count'))

count_visualizer(count_sets(data), 'pie', 'Set Count', os.path.join(STAT_DIR, 'Stat_count'))

stack_bar(data, 'Sets', range(2009, 2018), 'Year', os.path.join(STAT_DIR, 'sets_by_year.png'), split=';')

#count_visualizer(count_sets(data), 'tester', 'pie', 'Set Use')
"""

"""""
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