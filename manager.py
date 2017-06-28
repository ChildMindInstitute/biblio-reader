import os, pandas, sys, re, json


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

WORKING_PATH = dir(os.path.join(MAIN_DIR, 'working'))


# These are tuple sets that contain different regular expression terms Google Scholar matched and their categories,
# used to associate each publication with its attribute "Sets"

WEIGHTED_SETS = [('NKI', re.compile('(\Wnki\W*|nathan\s(s\. )?kline\sinstitute).*rockland|rockland\ssample')),
                 ('ADHD200', re.compile('adhd\W200')),
        ('CORR', re.compile('\scorr\s|(consortium\sfor\sreproducibility\sand\sreliability)')),
        ('ABIDE', re.compile('abide|(autism\sbrain\simaging\sdata\sexchange)'))]

UNWEIGHTED_SETS = [('FCP', re.compile('1,?000\sfunctional\sconnectomes?(\sproject)?|fcp|fcon[\s_]*1000'))]


# These are the papers that have connection to the main area of interest, found manually

CONTR_PAPERS = [1, 5, 74, 92, 653]


def get_file(file, dir):
    """This function retrieves the file based on the user-input filename and directory

    :param file: Name of the text file that contains the filename to be converted into a csv file for analysis
    :param dir: Directory name specified by WORKING_PATH
    :return: String containing filename to be converted to a csv file
    """
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
    """Converts the csv file specified in get_file into a pandas data structure

    :return: DATA -- contains a pandas data structure with Google Scholar results
    """
    if 'DATA' in globals():
        return
    global DATA_NAME, DATA
    DATA_NAME = get_file('DATA_NAME.txt', WORKING_PATH).read() + '.csv'
    DATA = pandas.read_csv(get_file(DATA_NAME, OUTPUT_PATH))
    return DATA


def get_paragraphs():
    """Retrieves text from the paragraphs.json file and converts into a dictionary

    :return: Dictionary with (key = integer index, value = paragraph text)
    """
    return {int(i): paragraph for i, paragraph in json.load(get_file('paragraphs.json', WORKING_PATH)).items()}


def get_bibs():
    """Use in affil_to_json.py to determine article affiliations and locate them geographically

    :return: Pandas data structure of article affiliations
    """
    return pandas.read_csv(get_file('parsed_bibs.csv', WORKING_PATH))


def update_data():
    print('data updated:', id(DATA))
    DATA.to_csv(path_or_buf=os.path.join(OUTPUT_PATH, DATA_NAME),
                index=False)


def get_journal_attrs():
    """Retrieves journal names and relevant attributes (Citescore, categories) from a json file located
    in inputs directory

    :return: A dictionary with (key = journal name, value = attributes)
    """
    return {key.lower(): value for key, value in json.load(get_file('journal_attrs.json', INPUT_PATH)).items()}