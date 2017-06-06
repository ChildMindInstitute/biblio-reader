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

def repair_affils(affiliations):
    aff_dict = dict()
    substitutions = [re.compile('\s\([^(]*\)'), re.compile('\s*(Electronic address:\s*)*\S+@\S+'),
                     re.compile('\A[^a-zA-Z]+'), re.compile('\s*[\.,]\Z'), re.compile('Email:.*')]
    for i, affs in affiliations.items():
        for aff in affs:
            for sub in substitutions:
                aff = re.sub(sub, '', aff)
            if re.search("([^./A-Z&'\s\d-])(?=[A-Z])", aff):
                aff = re.sub("([^./A-Z&'\s\d-])(?=[A-Z])", lambda x: x.group(0) + ' ' + x.group(1)[1:], aff)
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
    if os.path.exists('geo_affil_successes.txt'):
        with open('geo_affil_successes.txt', 'r') as a:
            geo_affils = a.readlines()
    else:
        geo_affils = []
    if os.path.exists('affiliations.json'):
        with open('affiliations.json', 'r') as js:
            geo_dict = {i: set(geo) for i, geo in json.load(js).items()}
    else:
        geo_dict = dict()
    for aff, ix in affiliations.items():
        if aff in geo_affils:
            continue
        request = geo_req(aff)
        try:
            geo_data = json.load(req.urlopen(request))
        except Exception as e:
            print(e, request)
            continue
        original_aff = str(aff)
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
        geo_affils.append(original_aff)
        print(aff, ': Success')
        location = geo_data['results'][0]['geometry']['location']
        latlong = ','.join([str(location['lat']), str(location['lng'])])
        if latlong not in geo_dict:
            geo_dict[latlong] = ix
        else:
            geo_dict[latlong].update(ix)
    with open('affiliations.json', 'w') as jf, open('geo_affil_successes.txt', 'w') as s:
        s.write('\n'.join(geo_affils))
        json.dump({latlong: [int(i) for i in ix] for latlong, ix in geo_dict.items()}, jf)

geo_lookup(repair_affils(affiliations))
