import random, collections, os, csv, ast
import pandas as pd

with open('../inputs/FCP_DATA.csv', 'r') as f, open('../outputs/paragraphs.txt', 'r') as p, \
        open('../outputs/unlinkables.csv', 'r') as u:
    data = pd.read_csv(f)
    paragraphs = ast.literal_eval(p.read())
    bad_data = pd.read_csv(u)


def checker_directory(directory):
    checks = {}
    for check in os.listdir(directory):
        full_path = '/'.join([directory, check])
        with open(full_path, 'r') as f:
            reader = list(csv.reader(f))
            for rows in reader[1:]:
                k = int(rows[0])
                v = (check, rows[1].replace(' and ', '').upper())
                if k not in checks:
                    checks[k] = [v]
                else:
                    checks[k].append(v)
    return checks


assign = list(data['i'])
checker = checker_directory('../inputs/Article_Checks')
double_checked = {key: check for key, check in checker.items() if len(check) > 1}
print(len([check for check in double_checked.values() if check[0][1] != check[1][1]]),
      *[check for check in double_checked.values() if check[0][1] != check[1][1]], sep='\n')
single_checked = [article for article in assign if article not in double_checked.keys()]


class Member(object):
    def __init__(self, name):
        self.name = name
        self.path = '../outputs/Assignments/' + name + '.txt'
        if os.path.exists(self.path):
            file = open(self.path)
            for i, line in enumerate(file):
                if i == 3:
                    self.articles = sorted([int(l) for l in str(line).strip().split(',')])
                    break
        else:
            self.articles = []

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in self.articles]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.iloc[key]['Title']) +
                        '\n' + str(data.iloc[key]['Authors']) + '\nPublication Type: ' +
                        data.iloc[key]['Journal Category'] + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraphs.items()) if key in self.articles]) + \
                '\n\n\n\n These numbers have no paragraphs:\n\n\n' + \
        '\n\n'.join(['ARTICLE NO ' + str(article) + ': ' + str(data.iloc[article]['Title']) +
                        '\n' + str(data.iloc[article]['Authors']) + '\n' + str(data.iloc[article]['URL']) for
                                   article in self.articles if article in list(bad_data['i'])])

    def assign(self, assign, length):
        while length > 0:
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                print(len(assign))
                length -= 1





class Assignment(object):
    def __init__(self, members):
        self.members = [Member(member) for member in members]

    def __getitem__(self, item):
        if item in range(0, len(self.members) - 1):
            return self.members[item]
        return None

    def assign(self, assignment, length=None):
        if length is None:
            length = len(assignment)
        length /= len(self.members)
        for member in self.members:
            member.assign(assignment, length)

    def write(self, new=None):
        for member in self.members:
            if new:
                with open(member.path.replace('.txt', '') + new + '.txt', 'w') as f:
                    f.write(str(member))
            else:
                with open(member.path, 'w') as f:
                    f.write(str(member))

assignment = Assignment(['anirudh_new', 'bonhwang_new', 'michael_new', 'helen_new', 'jon_new'])
