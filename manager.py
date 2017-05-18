import os, pandas
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = os.path.join(ROOT_DIR, 'inputs')
OUTPUT_PATH = os.path.join(ROOT_DIR, 'outputs')
DATA_PATH = os.path.join(INPUT_PATH, 'FCP_DATA.csv')
with open(DATA_PATH, 'r') as f:
    DATA = pandas.read_csv(f)


def update():
    DATA.to_csv(path_or_buf=DATA_PATH, index=False)