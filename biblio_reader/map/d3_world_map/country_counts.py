from collections import Counter
import os
import manager as mg
import pandas as pd

map_dir = mg.dir(os.path.join(mg.INPUT_PATH, 'map_tools'))
affiliations = list(mg.get_bibs()['affiliations'].dropna())

print(*['\n\n'.join(aff.split(';;')) for aff in affiliations], sep='\n\n')
states = [(state.strip().lower(), 'united states') for state in mg.get_file('states.txt', map_dir).readlines()]
countries = {tuple(line.strip().split(' | ')) for line in mg.get_file('countries.txt', map_dir).readlines()}
for country in countries:
    if len(country) == 1:
        print(country)
code_dict = dict(countries)
country_corrections = [tuple(corr.strip().split(' | ')) for corr in
                       mg.get_file('country_corrections.txt', map_dir).readlines()]
country_corrections += states

def getAffiliationCounts():
    counts = Counter()
    for aff in affiliations:
        aff = aff.strip().lower()
        aff_country = False
        for country in countries:
            country_name = country[0].lower()
            country_code = country[1].strip()
            if country_name in aff:
                counts[country_code] += 1
                aff_country = True
        corrected_countries = []
        for correction in country_corrections:
            if correction[0] in aff and correction[1] not in aff and \
                            correction[1] not in corrected_countries:
                counts[code_dict[correction[1].title()]] += 1
                corrected_countries.append(correction[1])
                aff_country = True
        if not aff_country:
            print("No country for", aff)
    return counts

def out(file):
    pd.DataFrame(list(getAffiliationCounts().items()), columns=['countries', 'counts']).to_csv(path_or_buf=file)

out(os.path.join(mg.ROOT_PATH, 'map_counts.csv'))
out('map_counts.csv')