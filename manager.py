import os, pandas, sys, ast


def dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    return dir


def write(file, directory, str):
    with open(os.path.join(dir(directory), file), 'w') as f:
        f.write(str)


MAIN_DIR = os.path.dirname(os.path.abspath(__file__))

INPUT_PATH = dir(os.path.join(MAIN_DIR, 'inputs'))

OUTPUT_PATH = dir(os.path.join(MAIN_DIR, 'outputs'))

ROOT_PATH = dir(os.path.join(MAIN_DIR, 'root'))


def get_data():
    if not os.path.isfile(os.path.join(ROOT_PATH, 'DATA.txt')):
        print('No data exists')
        sys.exit(1)
    with open(os.path.join(ROOT_PATH, 'DATA.txt'), 'r') as f:
        global DATA_NAME
        DATA_NAME = f.read() + '.csv'
    try:
        with open(os.path.join(OUTPUT_PATH, DATA_NAME), 'r') as d:
            global DATA
            DATA = pandas.read_csv(d)
            return DATA
    except IOError:
        print('Unable to open data file')
        sys.exit(1)

def get_paragraphs():
    if not os.path.isfile(os.path.join(ROOT_PATH, 'paragraphs.txt')):
        print('No paragraphs exist')
        sys.exit(1)
    try:
        with open(os.path.join(ROOT_PATH, 'paragraphs.txt'), 'r') as p:
            global PARAGRAPHS
            PARAGRAPHS = ast.literal_eval(p.read())
            return PARAGRAPHS
    except IOError:
        print('Unable to open paragraphs file')
        sys.exit(1)



def update():
    DATA.to_csv(path_or_buf=os.path.join(OUTPUT_PATH, DATA_NAME),
                index=False)

