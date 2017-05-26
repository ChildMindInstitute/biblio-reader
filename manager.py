import os, pandas, sys, ast, re


def dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir


def write(file, directory, st):
    if not isinstance(st, str):
        print('Not a string')
        sys.exit(1)
    with open(os.path.join(dir(directory), file), 'w') as f:
        f.write(st)


MAIN_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_PATH = dir(os.path.join(MAIN_DIR, 'inputs'))

OUTPUT_PATH = dir(os.path.join(MAIN_DIR, 'outputs'))

ROOT_PATH = dir(os.path.join(MAIN_DIR, 'root'))

TERMS = ['fcon_1000.projects.nitrc.org', 'Rockland Sample', '1000 Functional Connectomes',
       'International Neuroimaging Data-Sharing Initiative', 'Autism Brain Imaging Data Exchange', 'ADHD-200',
       'Consortium for Reproducibility and Reliability', 'FCP', 'ADHD 200', 'FCON 1000',
       'Functional Connectomes Project', 'www.nitrc.org/projects/fcon_1000', 'NITRC', 'ABIDE', 'Di Martino']


WEIGHTED_SETS = [('NKI', re.compile('(\WNKI\W)?Rockland\s(Sample)?')), ('ADHD200', re.compile('adhd\W200')),
        ('CORR', re.compile('\scorr\s|(consortium\sfor\sreproducibility\sand\sreliability)')),
        ('ABIDE', re.compile('abide|(autism\sbrain\simaging\sdata\sexchange)'))]

UNWEIGHTED_SETS = [('FCP', re.compile('1,?000\sfunctional\sconnectomes?(\sproject)?|fcp|fcon[\s_]1000')),
                   ('INDI', re.compile('\Windi\W|(international\sneuroimaging\sdata\Wsharing\sinitiative)'))]


def get_file(file, dir):
    file_name = file.replace('_', ' ').split('.')[0]
    if not os.path.isfile(os.path.join(dir, file)):
        print(file_name, 'does not exist')
        sys.exit(1)
    try:
        return open(os.path.join(dir, file))
    except IOError:
        print('Unable to open', file_name, 'file')
        sys.exit(1)

def get_data():
    global DATA_NAME, DATA
    DATA_NAME = get_file('DATA.txt', ROOT_PATH).read() + '.csv'
    DATA = pandas.read_csv(get_file(DATA_NAME, OUTPUT_PATH))
    return DATA

def get_paragraphs():
    return ast.literal_eval(get_file('paragraphs.txt', ROOT_PATH).read())


def get_bibs():
    return pandas.read_csv(get_file('parsed_bibs.csv', ROOT_PATH))


def update_data():
    DATA.to_csv(path_or_buf=os.path.join(OUTPUT_PATH, DATA_NAME),
                index=False)

def get_author_sets():
    df = pandas.read_csv(get_file('author_sets.csv', ROOT_PATH))
    return dict(zip(df['Author'], df['Data set']))
