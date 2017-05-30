import manager as mg
import json, re
from urllib import request as req
from unidecode import unidecode

bibs = mg.get_bibs().dropna(subset=['affiliations'])

API = 'AIzaSyCkBVxMHaiUrL1j-p_tc8fEFIbVxjjWqCk'
bibs['affiliations'] = bibs['affiliations'].apply(lambda aff: re.sub('\[?\d\]', ';', aff))
affiliations = {i: {aff.strip() for sublist in [affil.split(';') for affil in affiliation.split(';;')]
                    for aff in sublist} for i, affiliation in zip(bibs['i'], bibs['affiliations'])}

def repair_affils(affiliations):
    aff_dict = dict()
    substitutions = [re.compile('\s\([^(]*\)'), re.compile('\s*(Electronic address:\s*)*\S+@\S+'),
                     re.compile('\A[^a-zA-Z]+'), re.compile('\s*[\.,]\Z')]
    for i, affs in affiliations.items():
        for aff in affs:
            for sub in substitutions:
                aff = re.sub(sub, '', aff)
            if aff != '':
                if aff in aff_dict:
                    aff_dict[aff].add(i)
                else:
                    aff_dict[aff] = {i}
    return aff_dict

print(*repair_affils(affiliations).items(), sep='\n')

def geo_req(request):
    return 'https://maps.googleapis.com/maps/api/geocode/json?address=' + unidecode(request).replace(' ', '+') + \
            '&key=' + API


def geo_lookup(affiliations):
    geo_dict = dict()
    for aff, ix in affiliations.items():
        successes = 0
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
            continue
        successes += 1
        location = geo_data['results'][0]['geometry']['location']
        latlong = ','.join([str(location['lat']), str(location['lng'])])
        if latlong not in geo_dict:
            geo_dict[latlong] = [affils]
        else:
            geo_dict[latlong].append(affils)
        print('Successful finds for article no', affils, ':', successes)
    with open('affiliations.json', 'w') as jf:
        json.dump(geo_dict, jf)

#geo_lookup([])