import random, os, csv, ast, collections, math, manager
import pandas as pd

with open('../outputs/paragraphs.txt', 'r') as p, \
        open('../outputs/unlinkables.csv', 'r') as u:
    paragraphs = ast.literal_eval(p.read())
    bad_data = pd.read_csv(u)

data = manager.data
paragraphs = {key: value for key, value in paragraphs.items() if key not in [int(i) for i in bad_data['i']]}


def checker_directory(directory):
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
                    v = (check, rows[1].replace(' and ', '').upper())
                    if k not in checks:
                        checks[k] = [v]
                    elif len(checks[k]) < 2:
                        checks[k].append(v)
    return checks




double_checked = [(key, check) for key, check in
                  checker_directory('../inputs/validity_checks').items() if len(check) == 2]
print([i for i in range(0, 1560) if i not in [k[0] for k in double_checked]])
print(len(double_checked))
conflicts = [(key, check) for key, check in sorted(double_checked) if check[0][1] != check[1][1]]
conflicts_m = [(key, check) for key, check in conflicts if 'michael' in check[0][0] or 'michael' in check[1][0]]
conflicts_mj = [(key, check) for key, check in conflicts_m if 'jon' in check[0][0]]
conflicts_m = [(key, check) for key, check in conflicts_m if (key, check) not in conflicts_mj]
conflicts_no_m = [(key, check) for key, check in conflicts if not ('michael' in check[0][0] or 'michael' in check[1][0])]



class Member(object):
    def __init__(self, name, path):
        self.name = name
        self.path = '/'.join([path, name + '.txt'])
        if os.path.exists(self.path):
            file = open(self.path)
            for i, line in enumerate(file):
                if i == 3:
                    self.articles = sorted([int(l) for l in str(line).strip().split(',')])
                    break
        else:
            self.articles = []
        self.written = list(self.articles)

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in sorted(self.articles)]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(manager.iloc[key]['Title']) +
                        '\n' + str(manager.iloc[key]['Authors']) + '\nPublication Type: ' +
                        manager.iloc[key]['Journal Category'] + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraphs.items()) if key in self.articles]) + \
           '\n\n\n\nThese numbers have no paragraphs:\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(article) + ': ' + str(manager.iloc[article]['Title']) +
                        '\n' + str(manager.iloc[article]['Authors']) + '\n' + str(manager.iloc[article]['URL']) +
                        '\nPublication Type: ' + manager.iloc[article]['Journal Category'] for
                        article in self.articles if article in sorted(list(bad_data['i']))])

    def assign(self, assign, length):
        while length > 0 and not set(assign).issubset(set(self.articles)):
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1


class Assignment(object):
    def __init__(self, members=None, dir='../outputs/checker_assigns'):
        if members is None:
            self.members = [Member(member.replace('.txt', ''), dir) for member in os.listdir(dir) if
                            '.txt' in member]
        else:
            self.members = [Member(member, dir) for member in members]

    def __getitem__(self, item):
        if item in range(0, len(self.members)):
            return self.members[item]
        for member in self.members:
            if item == member.name:
                return member
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
                member.articles = [article for article in member.articles if article not in member.written]
                with open(member.path.replace('.txt', '') + new + '.txt', 'w') as f:
                    f.write(str(member))
            else:
                with open(member.path, 'w') as f:
                    f.write(str(member))

    def double(self):
        articles = []
        for member in self.members:
            articles += member.articles
        assignment = [item for item, count in collections.Counter(articles).items() if count == 1]
        member_length = math.ceil((len(articles) + len(assignment)) / len(self.members))
        revert = None
        while len(assignment) > 0:
            for member in self.members:
                if revert:
                     while revert > 0:
                         reverted = random.choice(member.articles)
                         if reverted not in assignment and reverted not in member.written:
                             assignment.append(reverted)
                             member.articles.remove(reverted)
                             revert -= 1
                     revert = None
                if len(member.articles) == member_length:
                    continue
                else:
                    member.assign(assignment, member_length - len(member.articles))
                    if len(member.articles) != member_length:
                        revert = len(assignment)

    def test(self, test_type=list([1, 2, 3])):
        if 1 in test_type:
            articles = []
            for member in self.members:
                articles += member.articles
            print('Duplicate items:',
                  len([item for item, count in collections.Counter(articles).items() if count == 2]))
        if 2 in test_type:
            for member in self.members:
                print(member.name.capitalize() + "'s Length:", len(member.articles), '| Duplicates:',
                      len([item for item, count in collections.Counter(member.articles).items() if count > 1]))
        if 3 in test_type:
            for member in self.members:
                print(member.name.capitalize() + "'s New Articles:",
                      len([article for article in member.articles if article not in member.written]), '| Old:',
                      len(member.written))


assignment = Assignment(members=['conflicts_m', 'conflicts_h', 'conflicts_j'])
assignment[0].articles = [key for key, check in conflicts_no_m]
assignment[1].articles = [key for key, check in conflicts_mj]
assignment[2].articles = [key for key, check in conflicts_m]
assignment.test()
assignment.write()