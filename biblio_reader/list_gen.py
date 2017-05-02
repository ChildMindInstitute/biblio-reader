import random, re, os
import pandas as pd
import collections

with open('../inputs/FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)


def find_paragraphs(txt_directory, terms):
    res = {}
    terms = list(map((lambda x: x.lower()), terms))
    for path, dirs, files in os.walk(txt_directory):
        for file in files:
            full_file = '/'.join([path, file])
            with open(full_file, 'r') as f:
                text = f.read()
            paragraphs = re.split(r'[ \t\r\f\v]*\n[ \t\r\f\v]*\n[ \t\r\f\v]*', text)
            paragraphs = [paragraph.lower() for paragraph in paragraphs if isinstance(paragraph, str)]
            key_paragraphs = []
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs += [paragraph for paragraph in paragraphs
                                   if re.search(re_term, paragraph) and paragraph not in key_paragraphs]
            for term in terms:
                re_term = term.replace(' ', '[\s]*')
                key_paragraphs = [str(re.sub(re_term, '@@@@' + term, paragraph)) for paragraph in key_paragraphs]
            res[int(file.replace('.txt', ''))] = key_paragraphs
    return res

fcp = ['fcon_1000.projects.nitrc.org', 'Rockland Sample', '1000 Functional Connectomes',
       'International Neuroimaging Data-Sharing Initiative', 'Autism Brain Imaging Data Exchange', 'ADHD-200',
       'Consortium for Reproducibility and Reliability', 'FCP', 'ADHD 200', 'FCON 1000',
       'Functional Connectomes Project', 'www.nitrc.org/projects/fcon_1000', 'NITRC']
paragraph_dict = find_paragraphs('../outputs/txts', fcp)
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
                    self.articles = [int(l) for l in str(line).strip().split(',')]
                    break
        else:
            self.articles = []

    def __str__(self):
        return self.name + '\n\n\n' + ','.join([str(article) for article in self.articles]) + '\n\n\n' + \
           '\n\n'.join(['ARTICLE NO ' + str(key) + ': ' + str(data.iloc[key]['Title']) +
                        '\n' + str(data.iloc[key]['Authors']) + '\n\n' + '\n\n'.join(paragraph)
                        for key, paragraph in sorted(paragraph_dict.items()) if key in self.articles])

    def assign(self, assign, length):
        while length > 0:
            assignment = random.choice(assign)
            if assignment not in self.articles:
                self.articles.append(assignment)
                assign.remove(assignment)
                length -= 1


class Assignment(object):
    def __init__(self, members):
        self.members = [Member(member) for member in members]

    def assign(self, assignment, length=None, duplicate=True):
        if length is None:
            length = len(assignment)
        length /= len(self.members)
        for member in self.members:
            if not duplicate:
                if len(member.articles) > length:
                    continue
            member.assign(assignment, length)

    def write(self):
        for member in self.members:
            with open(member.path, 'w') as f:
                f.write(str(member))

    def test(self):
        all_articles = []
        for member in self.members:
            all_articles += member.articles
        test_list = [item for item, count in collections.Counter(all_articles).items() if count > 1]
        print(len(test_list))
        print(test_list)