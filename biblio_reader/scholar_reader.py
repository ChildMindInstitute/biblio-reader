import pandas as pd
import matplotlib.pyplot as plt
import manager

def countstats(series, path=None, row_limit=None):
    stat = series.dropna().apply(lambda x: str(x).casefold()).value_counts()
    if row_limit:
        stat = stat.head(row_limit)
    if path:
        stat.to_csv(path=path, header=True)
    else:
        return stat


def citations_per_year(data, sort=False):
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


def count_item_year(data, column, item):
    return data.loc[data[column].apply(lambda x:str(x).casefold()) == item]['Year'].value_counts()


def year_sum_graph(data, column):
    df = data.groupby('Year').sum()[column]
    plt.figure()
    plt.ylabel(column)
    plt.xlabel('Year')
    df.plot.bar()
    plt.savefig('Citations_sum.png', bbox_inches='tight')




def categorize_journals(dict, keywords, category, type='Journal'):
    data_journals = data.dropna(subset=[type])
    for keyword in keywords:
        dict.update({row[1]['i']: category for row in data_journals.iterrows() if keyword in row[1][type].lower()})


other = ['arxiv', 'biorxiv', 'proceedings', 'advances in child development', 'using secondary datasets',
         '2 center for mind', 'under construct', 'bridging the gap before and after', 'autism imaging and devices']

journals = ['journal', 'plos', 'ieee', 'human brain mapping', 'frontiers', 'molecular', 'nature', 'neuro',
            'brain', 'biological psychiatry', 'autism', 'scientific', 'cortex', 'trends', 'jama', 'giga', 'bmc',
            'cell', 'chinese', 'phil. trans.', 'peerj', 'clinical trials', 'psychological', 'radiology', ]

journal_urls = ['nature', 'sciencedirect', 'wiley', 'springer', 'ncbi', 'journal']
other_urls = ['proceedings', 'ieeexplore', 'preprint', 'crcnetbase', 'patents']
t_urls = ['gradworks', 'academia', 'thesis', 'dissertation']

journal_categories = {}
categorize_journals(journal_categories, other, 'Other')
categorize_journals(journal_categories, journals, 'Journal')
categorize_journals(journal_categories, journal_urls, 'Journal', type='URL')
categorize_journals(journal_categories, other_urls, 'Other', type='URL')
categorize_journals(journal_categories, t_urls, 'Thesis/Dissertation', type='URL')
journal_categories.update({i: 'Unknown' for i in range(0, 1560) if i not in journal_categories})

#print(*sorted([(key, data.iloc[key]['Journal'], data.iloc[key]['URL']) for key, category in journal_categories.items() if category == 'Other']), sep='\n')
data['Journal Category'] = [journal for key, journal in sorted(journal_categories.items())]

print(*sorted(journal_categories.items()), sep='\n')
data.to_csv(path_or_buf='../inputs/FCP_DATA.csv', index=False)