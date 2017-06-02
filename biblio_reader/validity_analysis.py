import manager as mg
import os, sys, csv, collections
<<<<<<< HEAD
from biblio_reader.scholar_reader import categorize_journals
=======
from biblio_reader import scholar_reader
>>>>>>> 03143fbd02947324f278b2e2ca7f8438cc011736

checks = mg.dir(os.path.join(mg.INPUT_PATH, 'validity_checks'))
authors = mg.get_author_sets()

if len(os.listdir(checks)) == 0:
    print('No validity checks to analyze')
    sys.exit(1)

data = mg.get_data()
categories = mg.dir(os.path.join(mg.INPUT_PATH, 'journal_categories'))


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
<<<<<<< HEAD
    journal_types.update({key: typ for key, typ in categorize_journals(data, categories).items()
=======
    journal_types.update({key: typ for key, typ in scholar_reader.
                         categorize_journals(data, categories).items()
>>>>>>> 03143fbd02947324f278b2e2ca7f8438cc011736
                          if key not in journal_types})
    for key, type in journal_types.items():
        if type == 'Unknown':
            journal_types[key] = 'Journal'
    return {key: type for key, type in sorted(journal_types.items())}


def data_contributions_count(data, directory, author_associations):
    """
    If any manual investigator marks that a pub has some connection with the original source, the pub automatically
     gets marked as connected. All the papers with authors linked to that pub then get checked and added to the
     contributions count.

    :param data: The pandas dataframe
    :param directory: The directory of csv files of manual checks
    :param author_associations: A dictionary of authors and the terms that they are associated with
    :return: A list of papers that are considered part of the contributions count
    """
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
    print(len(original_contributers))
    contributing_authors = {author for i, authors in zip(data['i'], data['Authors']) for author in authors.split(' & ')
                            if i in contributing_papers}
    for row in data.dropna(subset=['Sets']).iterrows():
        row = row[1]
        authors = [author for author in row['Authors'].split(' & ') if author
                   in contributing_authors and author != 'others']
        sets = row['Sets'].split(';')
        i = row['i']
        if len(authors) != 0:
            for author in authors:
                print(author_associations[author], sets)
                if any(s in author_associations[author] for s in sets):
                    print(i, author, sets, author_associations[author])
                    contributing_papers.add(i)
    return contributing_papers
