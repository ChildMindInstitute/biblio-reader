import pandas as pd
import matplotlib.pyplot as plt
import manager as mg
import os, datetime, collections
import numpy as np
from biblio_reader import validity_analysis
STAT_DIR = mg.dir(os.path.join(mg.OUTPUT_PATH, 'stats'))
data = mg.get_data()


def count_visualizer(value_count, stat_type, name, row_limit=None):
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
    plt.savefig(os.path.join(STAT_DIR, name.lower().replace(' ', '_') + '.png'), bbox_inches='tight')
    plt.clf()


def citations_per_year(data, sort=False):
    """
    Calculates the number of citations per year for each pub
    :param data: The pandas dataframe
    :param sort: If true, re-indexes the dataframe based on citations per year
    """
    data['CPY'] = data['Citations'] / (datetime.datetime.now().year + 1 - data['Year'])
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)
        mg.update_data()


def stack_bar(data, column, stacker, stack_type, stat, split=None):
    """
    Almost the same as value counter, except each bar is stacked by a specific other column (such as finding out most
    popular journals by year, or term sets by usage, etc.) Examples are in the stats file
    :param data: The pandas dataframe
    :param column: The column of the dataframe to be value counted
    :param stacker: A list of distinct values that correspond to values in the stack_type
    :param stack_type: The column in the dataframe to be part of the stacks (such as year, etc.)
    :param split: If true, splits each row in the column by the splitter
    :return: A stacked bar graph
    """
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
    repaired_stacks = []
    for stack in stacks:
        empty_stacks = {key: 0 for key in max_stack.keys() if key not in stack}
        if len(empty_stacks) > 0:
            stack.update(empty_stacks)
            repaired_stacks.append({key: value for key, value in sorted(stack.items()) if key in max_stack})
        else:
            repaired_stacks.append(stack)
    stacks = repaired_stacks
    if stat == 'bar':
        last_stack = list(np.zeros(len(max_stack), dtype=np.int))
        for name, stack in zip(stacker, stacks):
            plt.bar(range(len(stack)), list(stack.values()), align='center', label=name,
                    bottom=last_stack)
            last_stack = [x + y for x, y in zip(last_stack, list(stack.values()))]
        plt.xticks(range(len(max_stack)), max_stack.keys())
    if stat == 'plot':
        plot_dict = dict()
        for stack in stacks:
            for key, value in stack.items():
                if key in plot_dict:
                    plot_dict[key].append(value)
                else:
                    plot_dict[key] = [value]
        for stack, plot in plot_dict.items():
            plt.plot(plot, label=stack)
        plt.xticks(range(len(stacker)), stacker)
    plt.legend()
    title = ' by '.join([column, stack_type])
    plt.title(title)
    plt.savefig(os.path.join(STAT_DIR, title.lower().replace(' ', '_') + '_' + stat + '.png'))

stack_bar(data[data['Data Use'] != 'Y'], 'Sets', range(2009, 2017), 'Year', 'plot', split=';')

def count_sets(data):
    """
    Takes the terms that Google Scholar matched for all the publications and counts how many of each there are
    :param data: The dataframe to be used
    :return: A value count series of terms
    """
    if 'Sets' not in data:
        raise IOError('No sets to work with')
    sets = [dset for dsets in data['Sets'].dropna() for dset in dsets.split(';')]
    return pd.Series(collections.Counter(sets), name='Sets')


def categorize_journals(data, categories):
    """
    DISCLAIMER: Must have an inputted list of keywords and corresponding categories to work with before using this
     function.

    Each publication that was found by Google Scholar is parsed through by keywords in Journal or URL columns to
    categorize each publication into specific categories

    :param data: The dataframe to use
    :param categories: A directory of text files named by category that contain keywords to categorize (see root)
    :return: A dictionary of publication indeces and their corresponding Journal Category
    """
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


def authors(data, link, split=None):
    """
    Links each author to every paper they are attributed to and to the specific column of that paper
    :param data: The pandas dataframe
    :param link: The column in the dataframe that will be linked in the result
    :param split: If not None, splits the link into its sublinks
    :return: A dictionary with author keys and a set of column links that they are linked to
    """
    author_links = {}
    data = data.dropna(subset=['Authors', link])
    for i, authors in zip(data['i'], data['Authors']):
        authors = [auth for auth in authors.split(' & ') if auth != 'others']
        for author in authors:
            linker = data.loc[i, link]
            if author in author_links:
                if split is not None:
                    author_links[author].update(set(linker.split(split)))
                else:
                    author_links[author].add(linker)
            else:
                if split is not None:
                    author_links[author] = set(linker.split(split))
                else:
                    author_links[author] = {linker}
    return author_links


def calculate_stats(data):
    """
    Prints the usage and journal category stats for the data set
    :param data: The pandas dataframe completed with data use and journal categories
    """
    if 'Data Use' not in data or 'Journal Category' not in data:
        print('Must first perform validity analysis')
        return
    contributions = validity_analysis.data_contributions_count\
        (data, mg.dir(os.path.join(mg.INPUT_PATH, 'validity_checks')), mg.get_author_sets())
    no_contributions = [i for i in range(0, len(data)) if i not in contributions]
    stats = []
    for usage in ['Y', 'S', 'N', 'I']:
        use_stats = []
        use_data = data[data['Data Use'] == usage]
        use_stats.append(len(use_data))
        use_stats.append(len(use_data[use_data['i'].isin(no_contributions)]))
        for type in ['Journal', 'Other', 'Thesis']:
            use_stats.append(len(use_data[use_data['Journal Category'] == type]))
        stats.append((usage, use_stats))
    for use_type, use in stats:
        for i, stat in enumerate(use):
            if i == 0:
                message = 'Number of publications'
            elif i == 1:
                message = 'Number of publications that did not contribute'
            elif i == 2:
                message = 'Number of journals'
            elif i == 3:
                message = 'Number of preprints, proceedings, books, etc,'
            else:
                message = 'Number of theses/dissertations'
            if use_type == 'Y':
                message += ' that used data:'
            elif use_type == 'N':
                message += ' that did not use data:'
            elif use_type == 'S':
                message += ' that only used scripts:'
            else:
                message += ' that were invalid:'
            print(message, stat)