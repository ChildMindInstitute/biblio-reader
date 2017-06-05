import manager as mg
import json, re, os
from urllib import request as req
from unidecode import unidecode

bibs = mg.get_bibs().dropna(subset=['affiliations'])
map_dir = mg.dir(os.path.join(mg.INPUT_PATH, 'map_tools'))
API = 'AIzaSyCkBVxMHaiUrL1j-p_tc8fEFIbVxjjWqCk'
bibs['affiliations'] = bibs['affiliations'].apply(lambda aff: re.sub('\[?\d\]', ';', aff))
affiliations = {i: {aff.strip() for sublist in [affil.split(';') for affil in affiliation.split(';;')]
                    for aff in sublist} for i, affiliation in zip(bibs['i'], bibs['affiliations'])}
valid_countries = [country.split(' | ')[0] for country in mg.get_file('countries.txt', map_dir).readlines()] + \
                  [state.strip() for state in mg.get_file('states.txt', map_dir).readlines()] + \
                  [corr.split(' | ')[0] for corr in mg.get_file('country_corrections.txt', map_dir).readlines()]
print(*valid_countries, sep='\n')

def repair_affils(affiliations):
    aff_dict = dict()
    substitutions = [re.compile('\s\([^(]*\)'), re.compile('\s*(Electronic address:\s*)*\S+@\S+'),
                     re.compile('\A[^a-zA-Z]+'), re.compile('\s*[\.,]\Z'), re.compile('Email:.*')]
    for i, affs in affiliations.items():
        for aff in affs:
            for sub in substitutions:
                aff = re.sub(sub, '', aff)
            valid = False
            for valids in valid_countries:
                if valids in aff:
                    valid = True
                    break
            if not valid:
                print(valid, aff)
            if aff != '':
                if aff in aff_dict:
                    aff_dict[aff].add(i)
                else:
                    aff_dict[aff] = {i}
    return aff_dict


def geo_req(request):
    return 'https://maps.googleapis.com/maps/api/geocode/json?address=' + unidecode(request).replace(' ', '+') + \
            '&key=' + API


def geo_lookup(affiliations):
    geo_dict = dict()
    for aff, ix in affiliations.items():
        request = geo_req(aff)
        try:
            geo_data = json.load(req.urlopen(request))
        except Exception as e:
            print(e, request)
            continue
        while geo_data['status'] == 'ZERO_RESULTS' and ',' in aff:
            aff = aff[aff.find(',') + 1:].strip()
            request = geo_req(aff)
            try:
                geo_data = json.load(req.urlopen(request))
            except Exception as e:
                print(e, request)
                break
        if len(geo_data['results']) == 0:
            print(aff, ': Failure')
            continue
        print(aff, ': Success')
        location = geo_data['results'][0]['geometry']['location']
        latlong = ','.join([str(location['lat']), str(location['lng'])])
        if latlong not in geo_dict:
            geo_dict[latlong] = ix
        else:
            geo_dict[latlong].update(ix)
    with open('affiliations.json', 'w') as jf, open('affiliations_fail.txt', 'w') as t:
        try:
            json.dump({latlong: [int(i) for i in ix] for latlong, ix in geo_dict.items()}, jf)
        except Exception as e:
            print(e)
            try:
                t.write(str(geo_dict))
            except:
                print(*geo_dict.items(), sep='\n')

#geo_lookup(repair_affils(affiliations))
affs = repair_affils(affiliations)

""""
for aff in affs.keys():
    if re.search(r'([^./A-Z\s\d(Mc)-])(?=[A-Z])', aff):
        print(aff)
        print(re.sub('([^./A-Z\s\d-])(?=[A-Z])', lambda x: x.group(0) + ' ' + x.group(1)[1:], aff), '\n\n\n')
"""