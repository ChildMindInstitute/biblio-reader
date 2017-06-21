import pandas as pd, matplotlib.pyplot as plt, manager as mg, os, datetime, collections, numpy as np, csv
from titlecase import titlecase
from scipy.stats import norm
STAT_DIR = mg.dir(os.path.join(mg.OUTPUT_PATH, 'stats'))
data = mg.get_data()


def count_visualizer(value_count, stat_type, name, row_limit=None, color=None):
    """
    Counts values of specific columns in dataframe
    :param value_count: A value counts series, dict, or LOT (see pandas value_count function)
    :param out: output file name
    :param stat_type: one of: bar, barh, pie, plot
    :param row_limit: Sets a limit to how many highest values should be counted
    :return: csv, bar, or pie file
    """
    value_count = {value: count for value, count in list(dict(value_count).items())[:row_limit]}
    plt.figure()
    if stat_type == 'bar':
        plt.bar(range(len(value_count)), list(value_count.values()), align='center', color=color)
        plt.xticks(range(len(value_count)), value_count.keys())
    elif stat_type == 'barh':
        plt.barh(range(len(value_count)), list(value_count.values()), align='center', tick_label=value_count.keys(),
                 color=color)
    elif stat_type == 'pie':
        plt.pie(list(value_count.values()), labels=value_count.keys(), autopct='%1.1f%%', shadow=True)
        plt.axis('equal')
    elif stat_type == 'plot':
        plt.plot(list(value_count.values()), color=color)
        plt.xticks([0, int(len(value_count) / 2), len(value_count)],
                   [list(value_count.keys())[0][1], list(value_count.keys())[int(len(value_count) / 2)][1],
                    list(value_count.keys())[len(value_count) - 1][1]])
        plt.fill_between(range(len(value_count)), list(value_count.values()))
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
    data['Citations Per Year'] = data['Citations'] / (datetime.datetime.now().year + 1 - data['Year'].apply(lambda x: (2000 + int(x[1:]))))
    if sort:
        data.sort_values('CPY', inplace=True, ascending=False)
        data.reset_index(drop=True, inplace=True)
        mg.update_data()


def journal_attrs(data, attr, count=False):
    data = data.dropna(subset=['Journal'])
    attrs = mg.get_journal_attrs()
    attrs = {journal: attrs[journal][attr] for journal in attrs}
    journals = sorted([(journal.lower(), attrs[journal.lower()]) for journal in data['Journal']
                       if journal.lower() in attrs], key=lambda attrib: attrib[1], reverse=True)
    if count:
        return collections.Counter(journals)
    else:
        return journals


def stacked_data(data, column, stack_type, stat, title=None, split=None, stacker_split=False):
    """
    Almost the same as value counter, except each type is stacked by a specific other column (such as finding out most
    popular journals by year, or term sets by usage, etc.) Examples are in the stats file
    :param data: The pandas dataframe
    :param column: The column of the dataframe to be value counted
    :param stack_type: The column in the dataframe to be part of the stacks (such as year, etc.)
    :param stat: One of: stacked, plot, cluster
    :param split: If true, splits each row in the column by the splitter
    :param stacker_split: If true, looks in each stack not for exactness but for inclusion
    :return: Either a stacked bar graph or a line plot, depending on the stat type
    """
    plt.figure()
    stacks = []
    stacker = list(data[stack_type].value_counts().index)[:15]
    for stack in stacker:
        if not split:
            if stacker_split:
                stacked = [value for value in data[data[stack_type].str.contains(stack).fillna(False)][column].dropna()]
            else:
                stacked = [value for value in data[data[stack_type] == stack][column].dropna()]
        else:
            if stacker_split:
                stacked = [value for values in
                           data[data[stack_type].str.contains(stack).fillna(False)][column].dropna()
                           for value in values.split(split)]
            else:
                stacked = [value for values in
                           data[data[stack_type] == stack][column].dropna()
                           for value in values.split(split)]
        stacks.append({titlecase(str(typ)): count for typ, count in sorted(collections.Counter(stacked).items())[:10]})
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
    if stat == 'stacked':
        last_stack = list(np.zeros(len(max_stack), dtype=np.int))
        for name, stack in zip(stacker, stacks):
            plt.bar(range(len(stack)), list(stack.values()), align='center', label=name,
                    bottom=last_stack)
            last_stack = [x + y for x, y in zip(last_stack, list(stack.values()))]
        plt.xticks(range(len(max_stack)), max_stack.keys())
    elif stat == 'plot':
        plot_dict = dict()
        for stack in stacks:
            for key, value in stack.items():
                if key in plot_dict:
                    plot_dict[key].append(value)
                else:
                    plot_dict[key] = [value]
        marker = 0
        for stack, plot in plot_dict.items():
            plt.plot(plot, label=stack, marker=marker)
            marker += 1
        plt.xticks(range(len(stacker)), stacker)
        column, stack_type = stack_type, column
    elif stat == 'cluster':
        curr_width = 0
        width = 0.75 / len(stacks)
        for name, stack in zip(stacker, stacks):
            plt.bar([x + curr_width for x in range(len(stack))], list(stack.values()), width, label=name)
            curr_width += width
        plt.xticks([x + (((len(stacks) / 2) - 0.5) * width) for x in range(len(max_stack))], max_stack.keys())
    else:
        raise IOError('Invalid stat type')
    plt.legend()
    if title is None:
        title = ' by '.join([column, stack_type])
    plt.title(title)
    plt.savefig(os.path.join(STAT_DIR, '_'.join([title.lower().replace(' ', '_'), stat]) + '.png'), bbox_inches='tight')



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
    :param categories: A directory of text files named by category that contain keywords to categorize (see working)
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
            res.update({i: cat_name for i, row in data_journals.iterrows() if keyword in row[type].lower()})
    res.update({i: 'Unknown' for i in range(len(data)) if i not in res})
    return {i: typ for i, typ in sorted(res.items())}


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


def impacts(data, sets, type, per_year=False):
    res = []
    for set in sets:
        if per_year:
            citations = sorted(data[data['Sets'].str.contains(set).fillna(False)]['Citations Per Year'], reverse=True)
        else:
            citations = sorted(data[data['Sets'].str.contains(set).fillna(False)]['Citations'], reverse=True)
        """
        plt.figure()
        plt.hist(citations, bins=max(citations), edgecolor='black')
        plt.xlabel('Total Number of Citations')
        plt.ylabel('Number of Publications')
        plt.title(set)
        plt.savefig(set + '.png')
        """

        if type == 'h':
            for i, citation in enumerate(citations):
                if citation >= i + 1:
                    h = citation
                else:
                    break
            res.append((set, h))
        elif type == 'mean':
            mean = sum(citations) / float(len(citations))
            res.append((set, mean))
        elif type == 'total':
            res.append((set, sum(citations)))
        elif type == 'sd':
            pass
            res.append((set, np.std(citations)))
        elif type == 'i10':
            i10 = len([cit for cit in citations if cit >= 10])
            res.append((set, i10))
        elif isinstance(type, int):
            res.append((set, np.percentile(citations, type)))
        else:
            raise IOError("Invalid type")
    return res


def calculate_stats(data):
    """
    Prints the usage, contributions, and journal category stats for the data set
    :param data: The pandas dataframe completed with data use and journal categories
    """
    if 'Data Use' not in data or 'Journal Category' not in data or 'Contributor' not in data:
        print('Must first perform validity analysis')
        return
    stats = []
    print(len(data[data['Contributor'] == 'Contributor']))
    for usage in ['Y', 'S', 'N', 'I']:
        use_stats = []
        use_data = data[data['Data Use'] == usage]
        use_stats.append(len(use_data))
        use_stats.append(len(use_data[use_data['Contributor'] == 'Not a Contributor']))
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


def data_contributions_count(data, update=False, original=False):
    """
    If any manual investigator marks that a pub has some connection with the original source, the pub automatically
     gets marked as connected

    :param data: The pandas dataframe
    :param directory: The directory of csv files of manual checks
    :param author_associations: A dictionary of authors and the terms that they are associated with
    :return: A list of papers that are considered part of the contributions count
    """
    contributing_papers = {1, 5, 74, 92, 653}
    if original:
        return contributing_papers
    author_associations = authors(data[data['i'].isin(contributing_papers)], 'Sets', split=';')
    for row in data.dropna(subset=['Sets']).iterrows():
        row = row[1]
        all_authors = [author for author in row['Authors'].split(' & ') if author
                   in author_associations and author != 'others']
        sets = row['Sets'].split(';')
        i = row['i']
        if len(all_authors) != 0:
            for author in all_authors:
                if any(s in author_associations[author] for s in sets):
                    contributing_papers.add(i)
    if update:
        data['Contributor'] = dict(sorted([(i, 'Contributor') for i in contributing_papers] +
                        [(i, 'Not a Contributor') for i in range(len(data)) if i not in contributing_papers])).values()
    return contributing_papers



#print(*sorted(data[data['Data Use'] == 'Y'][data['Sets'].str.contains(';').fillna(False)]['Citations']), sep='\n')
#print(*impacts(data[data['Data Use'] == 'Y'], count_sets(data).index, 'mean', per_year=True), sep='\n')
print(*impacts(data[data['Data Use'] == 'Y'], count_sets(data).index, 'sd'), sep='\n')