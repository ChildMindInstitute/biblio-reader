import random, os, csv, ast, collections
import pandas as pd

with open('../inputs/FCP_DATA.csv', 'r') as f, open('../outputs/paragraphs.txt', 'r') as p, \
        open('../outputs/unlinkables.csv', 'r') as u:
    data = pd.read_csv(f)
    paragraphs = ast.literal_eval(p.read())
    bad_data = pd.read_csv(u)

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
                k = int(rows[0])
                v = (check, rows[1].replace(' and ', '').upper())
                if k not in checks:
                    checks[k] = [v]
                else:
                    checks[k].append(v)
    return checks
double_checked = [(key, check) for key, check in checker_directory('../inputs/Article_Checks').items() if len(check) == 2]
conflicts = [(key, check) for key, check in sorted(double_checked) if check[0][1] != check[1][1]]

#print(len(double_checked))
#print(len(conflicts), *conflicts, sep='\n')

assign = list(data['i'])
assign_copy = list(assign)


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
        self.written = list(self.articles)

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in self.articles]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.iloc[key]['Title']) +
                        '\n' + str(data.iloc[key]['Authors']) + '\nPublication Type: ' +
                        data.iloc[key]['Journal Category'] + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraphs.items()) if key in self.articles]) + \
           '\n\n\n\nThese numbers have no paragraphs:\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(article) + ': ' + str(data.iloc[article]['Title']) +
                        '\n' + str(data.iloc[article]['Authors']) + '\n' + str(data.iloc[article]['URL']) +
                        '\nPublication Type: ' + data.iloc[article]['Journal Category'] for
                        article in self.articles if article in list(bad_data['i'])])

    def assign(self, assign, length):
        while length > 0 and not set(assign).issubset(set(self.articles)):
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1





class Assignment(object):
    def __init__(self, members):
        self.members = [Member(member) for member in members]


    def __getitem__(self, item):
        if item in range(0, len(self.members) - 1):
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

    def reassign(self):
        articles = []
        for member in self.members:
            articles += member.articles
        assignment = [item for item, count in collections.Counter(articles).items() if count == 1]
        member_length = (len(articles) + len(assignment)) / len(self.members)
        revert = None
        while len(assignment) > 0:
            for member in self.members:
                if revert:
                     while revert > 0:
                         reverted = random.choice(member.articles)
                         if reverted not in assignment:
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

    def test(self, test_type=None):
        if test_type == 1 or test_type is None:
            articles = []
            for member in self.members:
                articles += member.articles
            print('Duplicate items:',
                  len([item for item, count in collections.Counter(articles).items() if count == 2]))
        if test_type == 2 or test_type is None:
            for member in self.members:
                print(member.name.capitalize() + "'s Length:", len(member.articles), '| Duplicates:',
                      len([item for item, count in collections.Counter(member.articles).items() if count == 2]))





assignment = Assignment(['anirudh', 'bonhwang', 'michael', 'helen', 'jon'])

assignment.reassign()
assignment.reassign()
assignment.reassign()
assignment.test()

