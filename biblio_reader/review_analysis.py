import manager as mg, os, sys, csv, collections
data = mg.get_data()
from biblio_reader import scholar_reader
checks = mg.dir(os.path.join(mg.INPUT_PATH, 'article_review'))
categories = mg.dir(os.path.join(mg.INPUT_PATH, 'journal_categories'))
if len(os.listdir(checks)) == 0:
    print('No validity checks to analyze')
    sys.exit(1)


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
    global conflicts
    conflicts = [i for i, check in checks.items() if len(check) != 2]
    print('These numbers still have conflicts:', *conflicts, sep='\n')
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
    journal_types.update({key: typ for key, typ in scholar_reader.categorize_journals(data, categories).items()
                          if key not in journal_types})
    for key, type in journal_types.items():
        if type == 'Unknown':
            journal_types[key] = 'Journal'
    return {key: type for key, type in sorted(journal_types.items())}