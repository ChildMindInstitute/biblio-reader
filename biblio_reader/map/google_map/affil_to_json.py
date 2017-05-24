import manager as mg
import json, re
from urllib import request as req
from unidecode import unidecode

bibs = mg.get_bibs().dropna(subset=['affiliations'])

API = 'AIzaSyCkBVxMHaiUrL1j-p_tc8fEFIbVxjjWqCk'

affiliations = {i: {aff.strip() for sublist in [affil.split(';') for affil in affiliation.split(';;')]
                    for aff in sublist} for i, affiliation in zip(bibs['i'], bibs['affiliations'])}

def repair_affils(affiliations):
    for affils in affiliations:
        affs = affiliations[affils]
        repaired_affils = set()
        for aff in affs:
            aff_repaired = re.sub('\s\(.*\)', '', aff)
            aff_repaired = re.sub('\s*(Electronic address:\s*)*\S+@\S+', '', aff_repaired)
            aff_repaired = re.sub('\A[^a-zA-Z]+', '', aff_repaired)
            aff_repaired = re.sub('\s*[\.,]\Z', '', aff_repaired)
            if len(aff_repaired) != 0:
                repaired_affils.add(aff_repaired)
        affiliations[affils] = repaired_affils

repair_affils(affiliations)


def geo_req(request):
    return 'https://maps.googleapis.com/maps/api/geocode/json?address=' + unidecode(request).replace(' ', '+') + \
            '&key=' + API


def geo_lookup(affiliations):
    geo_dict = dict()
    for affils in affiliations:
        successes = 0
        affs = affiliations[affils]
        for aff in affs:
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

geo_lookup(affiliations)