import manager as mg
import json, re, os
from urllib import request as req
from unidecode import unidecode

bibs = mg.get_bibs().dropna(subset=['affiliations'])
map_dir = mg.dir(os.path.join(mg.INPUT_PATH, 'map_tools'))
API = 'AIzaSyCc3U_YDbluAh_Eja8Zc4e4PX04ndyDXgE'
bibs['affiliations'] = bibs['affiliations'].apply(lambda aff: re.sub('\[?\d\]', ';', aff))
affiliations = {i: {aff.strip() for sublist in [affil.split(';') for affil in affiliation.split(';;')]
                    for aff in sublist} for i, affiliation in zip(bibs['i'], bibs['affiliations'])}

def repair_affils(affiliations):
    aff_dict = dict()
    substitutions = [re.compile('\s\([^(]*\)'), re.compile('\s*(Electronic address:\s*)*\S+@\S+'),
                     re.compile('\A[^a-zA-Z]+'), re.compile('\s*[\.,]\Z'), re.compile('Email:.*'),
                     re.compile('\Aand\s*'), re.compile(',?\s+and\Z'), re.compile('tel:.*|.*affiliated.*|To whom.*')]
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
    if os.path.exists('affiliations.json'):
        with open('affiliations.json', 'r') as js:
            geo_dict = json.load(js)
        affils = {affil for latlong in geo_dict for affil in geo_dict[latlong]['affiliations']}
    else:
        geo_dict = dict()
        affils = set()
    for aff, ix in affiliations.items():
        if aff in affils:
            continue
        ix = {int(i) for i in ix}
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
        print(aff, ': Success')
        location = geo_data['results'][0]['geometry']['location']
        latlong = ','.join([str(location['lat']), str(location['lng'])])
        if latlong not in geo_dict:
            geo_dict[latlong] = {"papers": ix,
                                "affiliations": {original_aff},
                                "matched searches": {aff}}
        else:

            geo_dict[latlong]["papers"] = set(geo_dict[latlong]["papers"]).union(ix)
            geo_dict[latlong]["affiliations"] = set(geo_dict[latlong]["affiliations"]).union({original_aff})
            geo_dict[latlong]["matched searches"] = set(geo_dict[latlong]["matched searches"]).union({aff})
    geo_dict = {latlong: {"papers": list(geo_dict[latlong]['papers']), "affiliations":
        list(geo_dict[latlong]['affiliations']), "matched searches":
        list(geo_dict[latlong]['matched searches'])} for latlong in geo_dict}
    with open('affiliations.json', 'w') as jf:
        json.dump(geo_dict, jf)


