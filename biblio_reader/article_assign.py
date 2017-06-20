import random, os, collections, math, manager as mg, csv, matplotlib.pyplot as plt

checker_dir = mg.dir(os.path.join(mg.ROOT_PATH, 'checker_assigns'))
data = mg.get_data()
from biblio_reader import scholar_reader
data = data[data['i'].isin(scholar_reader.data_contributions_count(data))][data['Data Use'] == 'Y'][data['Sets'].str.contains('ABIDE').fillna(False) | data['Sets'].str.contains('ADHD200').fillna(False)]
paragraphs = mg.get_paragraphs()


class Member(object):
    def __init__(self, name, path):
        self.name = name
        self.path = '/'.join([path, name + '.txt'])
        if os.path.exists(self.path):
            file = open(self.path)
            if len(file.readlines()) < 4:
                self.articles = []
            else:
                for i, line in enumerate(file):
                    if i == 3:
                        self.articles = sorted([int(l) for l in str(line).strip().split(',')])
                        break
        else:
            self.articles = []
        self.written = list(self.articles)

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in sorted(self.articles)]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.loc[key, 'Title']) +
                        '\n' + str(data.loc[key, 'Authors']) + '\nPublication Type: ' +
                        data.loc[key, 'Journal Category'] + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraphs.items()) if key in self.articles])

    def assign(self, assign, length):
        while length > 0 and not set(assign).issubset(set(self.articles)):
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1


class Assignment(object):
    def __init__(self, members=None, dir=checker_dir):
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

"""
assignment = Assignment(members=['Michael', 'Anirudh', 'Jon', 'Bonhwang'], dir=mg.dir(os.path.join(checker_dir, 'contr_checks')))
assignment.assign(list(data['i']))
assignment.write()
"""

SIZES_DIR = os.path.join(mg.INPUT_PATH, 'sample_sizes')
sizes = {}
for size in os.listdir(SIZES_DIR):
    if '.csv' not in size:
        continue
    full_path = '/'.join([SIZES_DIR, size])
    with open(full_path, 'r') as f:
        reader = list(csv.reader(f))
        for rows in reader[1:]:
            sizes[int(rows[0])] = int(rows[1])

sizes = {article: size for article, size in sizes.items() if article in data['i']}
adhd200, abide = [], []
for article, size in sizes.items():
    if 'ABIDE' in data.loc[article, 'Sets'] and 'ADHD200' in data.loc[article, 'Sets']:
        abide.append(size / 2)
        adhd200.append(size/ 2)
    elif 'ABIDE' in data.loc[article, 'Sets']:
        abide.append(size)
    elif 'ADHD200' in data.loc[article, 'Sets']:
        adhd200.append(size)

print(abide, adhd200)

plt.figure()
plt.boxplot([abide, adhd200])
plt.xticks([1, 2], ['ABIDE', 'ADHD200'])
plt.title('Sample Sizes')
plt.savefig('tester.png')