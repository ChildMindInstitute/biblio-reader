import manager as mg
import os, sys, csv, collections
from biblio_reader import scholar_reader

checks = mg.dir(os.path.join(mg.INPUT_PATH, 'validity_checks'))
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


data['Journal Category'] = correct_types(checks, data).values()

data['Data Use'] = usage_directory(checks).values()

mg.update_data()

"""
double_checked_spec = [(key, check) for key, check in
                  specifier_directory('../inputs/validity_checks').items() if len(check) == 2]
cmi_authored = [key for key, check in double_checked_spec if 'Q' in check[0] or 'Q' in check[1]]

double_checked = [(key, check) for key, check in
                  checker_directory('../inputs/validity_checks').items() if len(check) == 2 and key not in cmi_authored]
conflicts = [(key, check) for key, check in double_checked if check[0] != check[1]]
usage = [key for key, check in double_checked if check[0] == 'Y']
no_usage = [key for key, check in double_checked if check[0] == 'N']
invalid = [key for key, check in double_checked if check[0] == 'I']
scripts = [key for key, check in double_checked if check[0] == 'S']
print(len(conflicts + usage + no_usage + invalid + scripts))
print(len(conflicts), len(usage), len(no_usage), len(invalid), len(scripts), len(double_checked))
"""