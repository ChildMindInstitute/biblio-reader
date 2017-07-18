#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Biblio_reader.working.vocab.controlled

This script reads and creates controlled vocabularies

Author:
    – Jon Clucas, 2017

Copyright ©2017, Child Mind Institute (http://childmind.org), Apache v2.0
    License
"""
import os
import sys
br_path = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir,
          os.pardir))
if br_path not in sys.path:
    sys.path.append(br_path)
from geopy.exc import GeocoderServiceError
from geopy.geocoders import GoogleV3 as Google
import json
import manager as mg
from geopy.geocoders import Nominatim
import biblio_reader.map.affil_to_json as atj
import pandas as pd
from termcolor import colored
import time
regappend = "(\s+|,+|$)"


def all_matched_searches(affiliations, de_facto_affiliations):
    """
    Function to create a dictionary for all matched geolocation search terms.

    Parameters
    ----------
    affiliations: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates; "country":
        dictionary of ISO strings representing countries

    de_facto_affiliations: set
        set of affiliation strings

    Returns
    -------
    affiliations: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates; "country":
        dictionary of ISO strings representing countries

    country_count: dictionary
        keys are Alpha-2 codes, values are counts
    """
    geolocator = Nominatim()
    backup_geolocator = Google("AIzaSyCc3U_YDbluAh_Eja8Zc4e4PX04ndyDXgE")
    iso_3166_1 = pd.read_csv(os.path.abspath(os.path.join(__file__, os.pardir,
                 "ISO_3166_1.csv")), na_filter=False)
    iso_3166_2_us = pd.read_csv(os.path.abspath(os.path.join(__file__,
                    os.pardir, "ISO_3166_2_US.csv")), na_filter=False)
    iso_dict = {**{country['Alpha-2 code']: [country[
               'English short name (upper/lower case)'], country[
               'Alpha-2 code'], country['Alpha-3 code']] for country in
               iso_3166_1.to_dict(orient='records')}, **{state['Code']: [
               state["Subdivision name"], state['Code'], state['Code']] for
               state in iso_3166_2_us.to_dict(orient='records')}, 'unknown': [
               'unknown'] * 3}
    countries = {**{country['Alpha-2 code']: country['Alpha-2 code'] for
                country in iso_3166_1.to_dict(orient='records')}, **{country[
                'Alpha-3 code']: country['Alpha-2 code'] for country in
                iso_3166_1.to_dict(orient='records')}, **{country[
                'English short name (upper/lower case)']: country[
                'Alpha-2 code'] for country in iso_3166_1.to_dict(orient=
                'records')}, **{state['Code']: state['Code'] for state in
                iso_3166_2_us.to_dict(orient='records')}, **{state[
                'Subdivision name']: state['Code'] for state in
                iso_3166_2_us.to_dict(orient='records')}, 'unknown': 'unknown',
                '?': 'unknown', 'Taiwan': 'TW', "PRC": "CN", "PR China": "CN",
                "UK": "GB", "United Kingdom": "GB"}
    us = {'US', 'USA', 'United States', 'U.S.A', "United States of America"}
    us_states = {state['Subdivision name']: state['Code'] for state in
                iso_3166_2_us.to_dict(orient='records')}
    for state in us_states:
        if state not in countries:
            countries[state] = us_states[state]
    country_count = {country: 0 for country in iso_dict}
    for k, v in affiliations.items():
        time.sleep(1)
        if "country" not in affiliations[k]:
            address_components = None
            while not address_components:
                time.sleep(1)
                try:
                    address_components = [x.strip() for x in
                                         geolocator.reverse(k, language=
                                         'en').address.split(',')]
                except GeocoderServiceError as g:
                    try:
                        address_components = list({com_g.strip() for com_g in [
                                             com_i for com_h in [com[0].split(
                                             ',') for com in
                                             backup_geolocator.reverse(k,
                                             language='en')] for com_i in com_h
                                             ]})
                    except:
                        print(colored(g, 'yellow'))
                        next
            if bool([u for u in us if u in address_components]):
                for state in us_states:
                    if state in address_components:
                        if state == "Georgia":
                            affiliations[k]["country"] = "US-GA"
                        else:
                            affiliations[k]["country"] = countries[state]
                        country_count[affiliations[k]["country"]] =           \
                                                                 country_count[
                                                                  affiliations[
                                                                    k][
                                                                    "country"]
                                                                    ] + 1
            else:
                for country in countries:
                    if country != 'United States of America' and country in   \
                       address_components:
                        affiliations[k]["country"] = countries[country]
                        country_count[affiliations[k]["country"]] =           \
                                                                 country_count[
                                                                  affiliations[
                                                                    k][
                                                                    "country"]
                                                                     ] + 1
            if "country" not in affiliations[k]:
                country = input(colored("{}\n{}? ".format(str(
                          address_components), str(affiliations[k][
                          "affiliations"])), 'magenta'))
                if len(country):
                    affiliations[k]["country"] = countries[country]
                    country_count[affiliations[k]["country"]] = country_count[
                                                                affiliations[
                                                                k]["country"]]\
                                                                + 1
        if "country" in affiliations[k]:
            print("{}: {}".format(iso_dict[affiliations[k]["country"]][0], str(
                  address_components)))
    save_heatmap_data(affiliations)
    return(affiliations, country_count)


def get_vocab(which_vocab):
    """
    Function to get the current cntrolled vocabulary for affiliations

    Parameters
    ----------
    which_vocab: string
        name of vocabulary to get

    Returns
    -------
    voc: dictionary
        dictionary with de facto spelling string keys and canonical spelling
        string values
    """
    path = os.path.join(mg.WORKING_PATH, 'vocab', ''.join([which_vocab, '.json'
           ]))
    if os.path.exists(path):
        with open(path, 'r') as js:
            return(json.load(js))
    else:
        return(dict())


def parse_affiliation(affiliation):
    """
    Function to parse affiliation strings

    Parameter
    ---------
    affiliation: string
        string to parse
    """
    aff_split = affiliation.split(', ')
    if len(aff_split) > 1:
        return(aff_split[0], ", ".join(aff_split[1:]))
    else:
        return(None, aff_split[0])


def parse_affiliations(path='affiliations.json'):
    """
    Function to parse existing affiliations data from file

    Parameter
    ---------
    path: string
        path to json, default='./affiliations.json'

    Returns
    -------
    geo_dict: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates

    affils: set
        set of de facto affiliations
    """
    geo_dict, affils = atj.get_affiliation_json(path)

    return geo_dict, affils


def save_heatmap_data(affiliations, outpath=os.path.abspath(os.path.join(
                      br_path, "biblio_reader", "map", "affiliations.csv"))):
    """
    Function to save the controlled heatmap counts

    Parameters
    ----------
    afiliations: dictionary
        dictionary with heatmap data

    outpath: string
        location to save data

    Returns
    -------
    None
    """
    country_count = dict()
    us_state_count = dict()
    for coord in affiliations:
        print(affiliations[coord]['country'])
        if affiliations[coord]['country'].startswith("US"):
            if affiliations[coord]['country'] not in us_state_count:
                us_state_count[affiliations[coord]['country']] = 0
            us_state_count[affiliations[coord]['country']] = us_state_count[
                           affiliations[coord]['country']] + len(affiliations[
                           coord]['papers'])
            affiliations[coord]['country'] = "US"
        if affiliations[coord]['country'] not in country_count:
            country_count[affiliations[coord]['country']] = 0
        country_count[affiliations[coord]['country']] = country_count[
                      affiliations[coord]['country']] + len(affiliations[coord
                      ]['papers'])

    country_count = pd.DataFrame(list(country_count.items()), columns=[
                      'country', 'count'])

    country_count[country_count['country'] != 'unknown'].to_csv(outpath, index=
                                                                False)

    us_state_count = pd.DataFrame(list(us_state_count.items()), columns=[
                      'country', 'count'])
    us_outpath = "".join([outpath.split("affiliations")[0],
                 "us_affiliations.csv"])
    us_state_count[us_state_count['country'] != 'unknown'].to_csv(us_outpath,
                                                                index=False)


def save_vocab(which_vocab, voc):
    """
    Function to save the current cntrolled vocabulary for affiliations

    Parameters
    ----------
    which_vocab: string
        name of vocabulary to get

    voc: dictionary
        dictionary with de facto spelling string keys and canonical spelling
        string values

    Returns
    -------
    None
    """
    path = os.path.join(mg.WORKING_PATH, 'vocab', ''.join([which_vocab, '.json'
           ]))

    cpath = "".join([path[:-5], ".csv"])

    v_path = os.path.dirname(path)
    if not os.path.exists(v_path):
        os.makedirs(v_path)

    with open(path, 'w') as js:
        json.dump(voc, js, sort_keys=True, indent=4)

    cvoc = pd.DataFrame()

    for i in voc:
        try:
            cvoc = cvoc.append(pd.DataFrame({which_vocab: i, **voc[i]}, index=[
                   len(cvoc)]))
        except:
            cvoc = cvoc.append(pd.DataFrame({which_vocab: i}, index=[len(cvoc)]
                   ))
    cvoc.to_csv(cpath, index=False)


def update_vocab(which_vocab):
    """
    Function to get the current cntrolled vocabulary for affiliations

    Parameters
    ----------
    which_vocab: string
        name of vocabulary to get

    Returns
    -------
    voc: dictionary
        dictionary with de facto spelling string keys and canonical spelling
        string values, updated

    count: dictionary
        dictionary of counts
    """
    up_fun = {'affiliation': (parse_affiliations, os.path.join(br_path,
             'biblio_reader', 'map', 'affiliations.json'))}
    geo_dict, new_v = up_fun[which_vocab][0](up_fun[which_vocab][1])
    voc, count = all_matched_searches(geo_dict, new_v)
    return(voc, count)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    which_vocab = 'affiliation'
    v = get_vocab(which_vocab)
    v, c = update_vocab(which_vocab)
    save_vocab(which_vocab, v)
    save_vocab('_'.join([which_vocab, 'count']), c)
