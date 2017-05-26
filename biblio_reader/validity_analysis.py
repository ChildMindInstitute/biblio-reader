import manager as mg
import os, sys, csv, collections
from biblio_reader import scholar_reader

checks = mg.dir(os.path.join(mg.INPUT_PATH, 'validity_checks'))
authors = mg.get_author_sets()

if len(os.listdir(checks)) == 0:
    print('No validity checks to analyze')
    sys.exit(1)

data = mg.get_data()
categories = mg.dir(os.path.join(mg.INPUT_PATH, 'journal_categories'))


def usage_directory(directory):
    checks = {}
    for check in os.listdir(directory):
        if '.csv' not in check:
            continue
        full_path = '/'.join([directory, check])
        with open(full_path, 'r') as f:
            reader = list(csv.reader(f))
            for rows in reader[1:]:
                if rows[1] != '':
                    k = int(rows[0])
                    v = rows[1].replace(' and ', '').upper()
                    if k not in checks:
                        checks[k] = [v]
                    else:
                        checks[k].append(v)
                    if len(checks[k]) > 2:
                        singles = [item for item, count in
                                   collections.Counter(checks[k]).items()
                                   if count == 1]
                        if len(singles) == 3:
                            continue
                        for single in singles:
                            checks[k].remove(single)
                        if len(checks[k]) == 3:
                            del checks[k][0]
    return {i: check[0] for i, check in sorted(checks.items()) if len(check) == 2}


def correct_types(directory, data):
    journal_types = {}
    for check in os.listdir(directory):
        if '.csv' not in check:
            continue
        full_path = '/'.join([directory, check])
        with open(full_path, 'r') as f:
            reader = list(csv.reader(f))
            for rows in reader[1:]:
                k = int(rows[0])
                v = rows[2].replace(' and ', '').upper()
                if 'T' in v:
                    journal_types[k] = 'Thesis'
                elif 'O' in v:
                    if k in journal_types:
                        if journal_types[k] != 'Thesis':
                            journal_types[k] = 'Other'
                    else:
                        journal_types[k] = 'Other'
    journal_types.update({key: typ for key, typ in scholar_reader.
                         categorize_journals(data, categories).items()
                          if key not in journal_types})
    for key, type in journal_types.items():
        if type == 'Unknown':
            journal_types[key] = 'Journal'
    return {key: type for key, type in sorted(journal_types.items())}


def data_contributions_count(data, directory, author_associations):
    data = data.dropna(subset=['Authors'])
    contributing_papers = set()
    for check in os.listdir(directory):
        if '.csv' not in check:
            continue
        full_path = '/'.join([directory, check])
        with open(full_path, 'r') as f:
            reader = list(csv.reader(f))
            for rows in reader[1:]:
                k = int(rows[0])
                v = rows[2].replace(' and ', '').upper()
                if 'Q' in v:
                    contributing_papers.add(k)
    global original_contributers
    original_contributers = list(contributing_papers)
    contributing_authors = {author for i, authors in zip(data['i'], data['Authors']) for author in authors.split(' & ')
                            if i in contributing_papers}
    for row in data.iterrows():
        row = row[1]
        authors = [author for author in row['Authors'].split(' & ') if author in contributing_authors]
        sets = row['Sets']
        i = row['i']
        if len(authors) != 0:
            for author in authors:
                if any(s in author_associations[author] for s in sets):
                    contributing_papers.add(i)
    return contributing_papers


def checker_stat(type, checks, intersection=None):
    stat = [key for key, check in checks.items() if check[0] == type]
    if not intersection:
        return len(stat)
    else:
        return len([key for key in stat if key in intersection])

def calculate_stats():
    use_checks = usage_directory(checks)
    type_checks = correct_types(checks, data)
    contributions = data_contributions_count(data, checks, authors)
    no_contributions = [i for i in range(0, len(data)) if i not in contributions]
    stats = []
    for usage in ['Y', 'S', 'N', 'I']:
        use_stats = []
        use_stats.append(checker_stat(usage, use_checks))
        use_stats.append(checker_stat(usage, use_checks, intersection=no_contributions))
        for type in ['Journal', 'Other', 'Thesis']:
            use_stats.append(checker_stat(usage, use_checks,
                                          intersection=[key for key, check in type_checks.items() if check == type]))
        stats.append((usage, use_stats))
    for use_type, use in stats:
        for i, stat in enumerate(use):
            if i == 0:
                message = 'Number of articles'
            elif i == 1:
                message = 'Number of articles that did not contribute'
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

"""
def calculate_stats():
    pass
    usage = usage_directory(checks)
    types = correct_types(checks, data)
    contributions = data_contributions_count(data, checks, authors)
    no_contributions = [i for i in range(0, len(data)) if i not in contributions]
    print('Number of articles that used data:', checker_stat('Y', usage), '\n',
          'Number of articles that did not use data:', checker_stat('N', usage), '\n',
          'Number of articles that only used scripts:', checker_stat('S', usage), '\n',
          'Number of articles that were invalid:', checker_stat('I', usage))
    print('Number of journals:', checker_stat('Journal', types), '\n',
          'Number of theses and dissertations:', checker_stat('Thesis', types), '\n',
          'Number of preprints, proceedings, books, etc.:', checker_stat('Other', types))
    print('Number of articles that used data and did not contribute:', checker_stat('Y', usage, intersection=no_contributions))
"""

calculate_stats()
"""
conflicts = [(key, check) for key, check in double_checked if check[0] != check[1]]
usage = [key for key, check in double_checked if check[0] == 'Y']
no_usage = [key for key, check in double_checked if check[0] == 'N']
invalid = [key for key, check in double_checked if check[0] == 'I']
scripts = [key for key, check in double_checked if check[0] == 'S']
print(len(conflicts + usage + no_usage + invalid + scripts))
print(len(conflicts), len(usage), len(no_usage), len(invalid), len(scripts), len(double_checked))
"""