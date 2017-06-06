import manager as mg
import os, sys, csv, collections

checks = mg.dir(os.path.join(mg.INPUT_PATH, 'validity_checks'))

if len(os.listdir(checks)) == 0:
    print('No validity checks to analyze')
    sys.exit(1)

data = mg.get_data()
print(id(data))
categories = mg.dir(os.path.join(mg.INPUT_PATH, 'journal_categories'))


def categorize_journals(data, categories):
    from biblio_reader import scholar_reader
    return scholar_reader.categorize_journals(data, categories)


def author_links(data, link, split=None):
    from biblio_reader import scholar_reader
    return scholar_reader.authors(data, link, split=split)

def usage_directory(directory):
    """
    After manually looking at each publication and marking correctly, takes the marks from the csv directory
     and makes sure each are double checked for accuracy.
    :param directory: The directory of csv files of manual checks
    :return: A dictionary of the publication numbers and their validity
    """
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
    print('These numbers still have conflicts:', [i for i, check in checks.items() if len(check) != 2])
    return {i: check[0] for i, check in sorted(checks.items()) if len(check) == 2}


def correct_types(directory, data):
    """
    After manually going through each pub and figuring out what type of pub it is, takes the manual and automatic parsing
    And creates a dictionary of publications and their publication type

    1. If any manual investigator marks a pub as Thesis or Dissertation, the pub is automatically categorized as a Thesis
    2. If any manual investigator marks a pub as Other and no one marks Thesis, it is categorized as Other
    3. All other pubs are categorized by the automatic parser.

    :param directory: The directory of csv files of manual checks
    :param data: The pandas dataframe
    :return: A dictionary of publications and their publication type
    """
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
    journal_types.update({key: typ for key, typ in categorize_journals(data, categories).items()
                          if key not in journal_types})
    for key, type in journal_types.items():
        if type == 'Unknown':
            journal_types[key] = 'Journal'
    return {key: type for key, type in sorted(journal_types.items())}

def data_contributions_count(data, directory, update=False):
    """
    If any manual investigator marks that a pub has some connection with the original source, the pub automatically
     gets marked as connected. All the papers with authors linked to that pub then get checked and added to the
     contributions count.

    :param data: The pandas dataframe
    :param directory: The directory of csv files of manual checks
    :param author_associations: A dictionary of authors and the terms that they are associated with
    :return: A list of papers that are considered part of the contributions count
    """
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
    author_associations = author_links(data[data['i'].isin(contributing_papers)], 'Sets', split=';')
    for row in data.dropna(subset=['Sets']).iterrows():
        row = row[1]
        authors = [author for author in row['Authors'].split(' & ') if author
                   in author_associations and author != 'others']
        sets = row['Sets'].split(';')
        i = row['i']
        if len(authors) != 0:
            for author in authors:
                if any(s in author_associations[author] for s in sets):
                    contributing_papers.add(i)
    if update:
        data['Contributor'] = dict(sorted([(i, 'Contributor') for i in contributing_papers] +
                        [(i, 'Not a Contributor') for i in range(len(data)) if i not in contributing_papers])).values()
    print(len(contributing_papers))
    return contributing_papers
#data_contributions_count(data, checks, update=True)
del data['Contributor']
mg.update_data()