#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates a JSON file with affiliation information for mapmaking

Authors:
    – Michael Fleischmann, 2017
    – Jon Clucas, 2017

Copyright ©2017, Child Mind Institute (http://childmind.org), Apache v2.0
    License
"""
import os
import sys
br_path = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(
          __file__))))
if br_path not in sys.path:
    sys.path.append(br_path)
import manager as mg
import json, re
from urllib import request as req
from unidecode import unidecode


def main():
    bibs = mg.get_bibs().dropna(subset=['affiliations'])
    map_dir = mg.dir(os.path.join(mg.INPUT_PATH, 'map_tools'))
    API = 'AIzaSyCc3U_YDbluAh_Eja8Zc4e4PX04ndyDXgE'
    bibs['affiliations'] = bibs['affiliations'].apply(lambda aff: re.sub(
                           '\[?\d\]', ';', aff))
    affiliations = {i: {aff.strip() for sublist in [affil.split(';') for affil
                   in affiliation.split(';;')] for aff in sublist} for i,
                   affiliation in zip(bibs['i'], bibs['affiliations'])}
    parse_affiliations()


def get_affiliation_json(path='affiliations.json'):
    """
    Function to read existing affiliation data from file

    Parameter
    ---------
    path: string
        path to json, default='./affiliations.json'

    Returns
    -------
    2-tuple:

    geo_dict: dictionary
        dictionary from json file

    affils: set
        set of affiliations (strings) from geo_dict
    """
    if os.path.exists(path):
        with open(path, 'r') as js:
            geo_dict = json.load(js)
        affils = {affil for latlong in geo_dict for affil in geo_dict[latlong][
                 'affiliations']}
    else:
        geo_dict = {}
        affils = set()
    return(geo_dict, affils)


def repair_affils(affiliations):
    """
    Function to remove special characters and electronic addresses (i.e.,
    emails and phone numbers) from affiliation strings.

    Parameter
    ---------
    affiliations: dictionary
        dictionary of key, affiliation (string) pairs

    Returns
    -------
    aff_dict: dictionary
        dictionary of key, affiliation (string) pairs
    """
    aff_dict = dict()
    substitutions = [re.compile('\s\([^(]*\)'), re.compile(
                    '\s*(Electronic address:\s*)*\S+@\S+'), re.compile(
                    '\A[^a-zA-Z]+'), re.compile('\s*[\.,]\Z'), re.compile(
                    'Email:.*'), re.compile('\Aand\s*'), re.compile(
                    ',?\s+and\Z'), re.compile('tel:.*|.*affiliated.*|To whom.*'
                    )]
    for i, affs in affiliations.items():
        for aff in affs:
            for sub in substitutions:
                aff = re.sub(sub, '', aff)
            if re.search("([^./A-Z&'\s\d-])(?=[A-Z])", aff):
                aff = re.sub("([^./A-Z&'\s\d-])(?=[A-Z])", lambda x: x.group(0)
                      + ' ' + x.group(1)[1:], aff)
            if aff != '':
                if aff in aff_dict:
                    aff_dict[aff].add(i)
                else:
                    aff_dict[aff] = {i}
    return aff_dict


def geo_req(request):
    """
    Function to insert a geograpdic search term into Google Maps' REST API.

    Parameter
    ---------
    request: string
        search term

    Returns
    -------
    url: string
        RESTful URL
    """
    return 'https://maps.googleapis.com/maps/api/geocode/json?address=' +     \
           unidecode(request).replace(' ', '+') + '&key=' + API


def geo_lookup(affiliations):
    """
    Function to lookup missing geodata and append said data into appropriate
    JSON file.

    Parameter
    ---------
    affilitions: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates

    Output
    ------
    affiliations.json: JSON file
        updated version of file if extant, otherwise new file; keys are
        latitude, longitude pairs as strings; values are dictionaries of
        "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates


    Returns
    -------
    None
    """
    geo_dict, affils = get_affiliation_json()
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

            geo_dict[latlong]["papers"] = set(geo_dict[latlong]["papers"]
                                          ).union(ix)
            geo_dict[latlong]["affiliations"] = set(geo_dict[latlong][
                                                "affiliations"]).union({
                                                original_aff})
            geo_dict[latlong]["matched searches"] = set(geo_dict[latlong][
                                                    "matched searches"]).union(
                                                    {aff})
    geo_dict = {latlong: {"papers": list(geo_dict[latlong]['papers']),
               "affiliations": list(geo_dict[latlong]['affiliations']),
               "matched searches": list(geo_dict[latlong]['matched searches'])}
               for latlong in geo_dict}
    with open('affiliations.json', 'w') as jf:
        json.dump(geo_dict, jf)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
