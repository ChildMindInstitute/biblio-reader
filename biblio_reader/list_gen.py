import random, re, os
import pandas as pd

with open('../inputs/FCP_DATA.csv', 'r') as f:
    data = pd.read_csv(f)


fcp = ['fcon_1000.projects.nitrc.org', 'Rockland Sample', '1000 Functional Connectomes',
       'International Neuroimaging Data-Sharing Initiative', 'Autism Brain Imaging Data Exchange', 'ADHD-200',
       'Consortium for Reproducibility and Reliability', 'FCP', 'ADHD 200', 'FCON 1000',
       'Functional Connectomes Project', 'www.nitrc.org/projects/fcon_1000', 'NITRC']
jon, michael, bonhwang, helen, anirudh = [], [], [], [], []
assign = list(data['i'])
url_dict = dict(zip(data['i'], data['URL']))
valid_urls = {key: url for key, url in url_dict.items() if not isinstance(url, float)}
no_paragraphs = [1011, 1071, 1075, 1077, 1094, 1098, 1100, 1171, 1199, 1243, 1260, 1272, 1284, 1304, 1343, 1368, 1433,
                 1448, 1450, 1464, 1517, 1557, 201, 24, 245, 314, 332, 357, 363, 473, 496, 555, 654, 73, 76, 774, 996]
pdfs = [int(pdf.replace('.pdf', '')) for pdf in os.listdir('../inputs/pdfs') if not pdf.startswith('.')]
no_pdfs = [assignment for assignment in assign if assignment not in pdfs]
jon += no_paragraphs + no_pdfs
assign = [a for a in assign if a not in jon]
assign_copy = list(assign)


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
                term = term.replace(' ', '[\s]*')
                key_paragraphs += [paragraph for paragraph in paragraphs
                                   if re.search(term, paragraph) and paragraph not in key_paragraphs]
                key_paragraphs = [paragraph.replace(term, '@@@@' + term) for paragraph in key_paragraphs]
            res[int(file.replace('.txt', ''))] = key_paragraphs
    return res


paragraph_dict = find_paragraphs('../outputs/txts', fcp)


def assign_articles(member, assign, length):
    while length > 0:
        assignment = random.choice(assign)
        if assignment not in member:
            member.append(assignment)
            assign.remove(assignment)
            length -= 1


def assign_paragraphs(member):
    return '\n\n'.join(['ARTICLE NO ' + str(key) + '\n\n' + '\n\n'.join(paragraph)
                 for key, paragraph in sorted(paragraph_dict.items()) if key in member])


def assignment(write=True):
    assign_articles(michael, assign, 312)
    assign_articles(jon, assign, 312 - len(jon))
    assign_articles(helen, assign, 312)
    assign_articles(bonhwang, assign, 312)
    assign_articles(anirudh, assign, 312)
    print([i for i in range(0, 1559) if i not in michael + jon + helen + bonhwang + anirudh],
          len(michael + jon + helen + bonhwang + anirudh))
    assign_articles(michael, assign_copy, 80)
    assign_articles(jon, assign_copy, 80)
    assign_articles(helen, assign_copy, 80)
    assign_articles(bonhwang, assign_copy, 80)
    assign_articles(anirudh, assign_copy, 80)
    if write:
        with open('../outputs/Assignments/michael.txt', 'w') as m:
            m.write(assign_paragraphs(michael))
        with open('../outputs/Assignments/helen.txt', 'w') as h:
            h.write(assign_paragraphs(helen))
        with open('../outputs/Assignments/jon.txt', 'w') as j:
            j.write(assign_paragraphs(jon))
        with open('../outputs/Assignments/bonhwang.txt', 'w') as b:
            b.write(assign_paragraphs(bonhwang))
        with open('../outputs/Assignments/anirudh.txt', 'w') as a:
            a.write(assign_paragraphs(anirudh))

