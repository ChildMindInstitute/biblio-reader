import re, os, json
import manager as mg
from biblio_reader.text_tools import convertToText


def find_paragraphs(txt_directory, terms, outfile=None):
    """
    Finds the exact paragraphs in which the terms entered into Google Scholar were found.

    :param txt_directory: The directory of txt files converted from pdfs found from pdf_utils
    :param terms: A list of tuples containing the terms that were searched and a regex pattern that can be used to find
                    Each term
    :param outfile: If not None, outputs a json file with each distinct pub index number and its corresponding paragraphs
    :return: The dictionary with pub index numbers and corresponding paragraphs
    """
    res = {}
    # terms = list(map((lambda x: x.lower()), terms))
    for file in os.listdir(txt_directory):
        full_file = os.path.join(txt_directory, file)
        with open(full_file, 'r') as f:
            text = f.read()
        paragraphs = re.split(r'[ \t\r\f\v]*\n[ \t\r\f\v]*\n[ \t\r\f\v]*', text)
        paragraphs = [paragraph.lower() for paragraph in paragraphs if isinstance(paragraph, str)]
        key_paragraphs = []
        for term, pattern in terms:

            key_paragraphs += {paragraph for paragraph in paragraphs if re.search(pattern, paragraph)}
        for term, pattern in terms:
            key_paragraphs = [str(re.sub(pattern, '@@@@' + term, paragraph)) for paragraph in key_paragraphs]
        res[int(file.replace('.txt', ''))] = key_paragraphs
    if outfile:
        with open(outfile + '.json', 'w') as f:
            json.dump(res, f)
    return res


def assoc_sets(data, dir, sets, less_weighted_sets=None):
    """
    Sets are ways to categorize the search terms used in the initial Google Scholar search. Each full text is parsed to
    try and find which search terms Google Scholar matched, and those terms are matched to each article reference number
    :param dir: The directory of txt files converted from pdfs found from pdf_utils
    :param sets: A list of tuples that contain categories of terms and the regular expressions that Google Scholar might
     have found for each category. These terms which will be looked at first and associated first
    :param less_weighted_sets: An optional less weighted set of terms that will be associated if none come up for the
                                first set of terms
    :return: A dict of index numbers and their corresponding terms
    """
    res = {}
    for file in os.listdir(dir):
        associations = []
        full_file = os.path.join(dir, file)
        with open(full_file, 'r') as f:
            text = f.read().lower()
        for group, pattern in sets:
            if re.search(pattern, text) is not None:
                associations.append(group)
        if less_weighted_sets is not None and len(associations) == 0:
            for group, pattern in less_weighted_sets:
                if re.search(pattern, text) is not None:
                    associations.append(group)
                    break
        res[int(file.replace('.txt', ''))] = associations
    res.update({i: [] for i in range(len(data)) if i not in res})
    res = {i: ';'.join(sets) for i, sets in sorted(res.items())}
    return res


def main():
    if 'WEIGHTED_SETS' not in dir(mg) and 'UNWEIGHTED_SETS' not in dir(mg):
        print('Define WEIGHTED_SETS and UNWEIGHTED_SETS in manager.py before using this module')
        return
    TXT_DIR = mg.dir(os.path.join(mg.WORKING_PATH, 'txts'))
    PDF_DIR = mg.dir(os.path.join(mg.INPUT_PATH, 'pdfs'))
    data = mg.get_data()
    convertToText.walkAndText(PDF_DIR, TXT_DIR)
    find_paragraphs(TXT_DIR, mg.WEIGHTED_SETS + mg.UNWEIGHTED_SETS, outfile=os.path.join(mg.WORKING_PATH, 'paragraphs'))
    sets = assoc_sets(data, TXT_DIR, mg.WEIGHTED_SETS, less_weighted_sets=mg.UNWEIGHTED_SETS)
    data['Sets'] = sets.values()
    mg.update_data()

