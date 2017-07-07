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
import csv
import json
import manager as mg
import biblio_reader.map.affil_to_json as atj
import pandas as pd


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
        dictionary of strings representing countries as spelled and canonically

    de_facto_affiliations: set
        set of affiliation strings

    Returns
    -------
    affiliations: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates; "country":
        dictionary of strings representing countries as spelled and canonically
    """
    ll = dict()
    iso_3166_1_en_short = list()
    with open("ISO_3166_1_English_short_names.csv", "r") as iso:
        iso_reader = csv.reader(iso)
        for line in iso_reader:
            iso_3166_1_en_short += line
    countries = dict()
    for s_name in iso_3166_1_en_short:
        countries[s_name] = s_name
    for coord, v in affiliations.items():
        for term in v["matched searches"]:
            ll[term] = coord
    latlong = dict()
    for affil in de_facto_affiliations:
        if affil in ll:
            latlong[affil] = {"match": affil, "coords": ll[affil]}
        else:
            aff_l, aff_r = parse_affiliation(affil)
            if aff_r in ll:
                latlong[affil] = {"match": aff_r, "coords": ll[aff_r]}
            else:
                while(aff_r not in ll and len(aff_r.split(", ")) > 1):
                    aff_l, aff_r = parse_affiliation(aff_r)
                if(aff_r in ll):
                    latlong[affil] = {"match": aff_r, "coords": ll[aff_r]}
                else:
                    latlong[affil] = None
    for affil in latlong:
        matched_country = None
        canon_match = None
        if affil and "country" not in latlong[affil]:
            while not matched_country:
                matched_country = country_prompt(countries, affil)
            if matched_country in countries:
                latlong[affil]["country"] = {'de facto': matched_country,
                                             'canonical': countries[
                                              matched_country]}
            else:
                if(len(countries)):
                    print("\n")
                    for cc in sorted(list({cc for c, cc in countries.items()})
                              ):
                        print(cc, end="; ")
                    print("\n")
                    print("If that country is the same as one of the above,"
                          " please enter the match. Otherwise, just press "
                          "`enter`")
                    canon_match = input("(please type the full name of the "
                                  "match, matching case): ")
                canon_match = matched_country if not canon_match else         \
                              canon_match
                if matched_country not in affil:
                    matched_country = input("How's that spelled here? ")
                if not matched_country or len(matched_country) == 0:
                    matched_country = canon_match
                latlong[affil]["country"] = {'de facto': matched_country,
                                             'canonical': canon_match}
            if canon_match:
                countries[matched_country] = canon_match
                if canon_match not in countries:
                    countries[canon_match] = canon_match
    for coord, v in affiliations.items():
        for term in v["matched searches"]:
            affiliations[coord]["country"] = latlong[term]["country"]
    return(affiliations)


def country_prompt(countries, affil):
    """
    Function to prompt user if affiliation is in a given country, and, if not,
    to identify the affiliation.

    Parameters
    ----------
    countries: dictionary
        countries vocabulary

    affil: string
        affiliation

    Returns
    -------
    country: string or None
        confirmed country
    """
    doublecheck = ["India", "CA", "MA", "Georgia", "US", "Mexico"]
    for country in countries:
        if country and country in affil:
            if country not in doublecheck:
                return(country)
            match_country = "maybe"
            while(len(match_country) > 0 and match_country.upper()
                  not in ["Y", "N"]):
                match_country = input(''.join(["Is ", affil, " in the country",
                                " \"", countries[country], "\"? (Y / N) : "]))
            if len(match_country) == 0 or match_country.upper() == "Y":
                return(country)
    match_country = input(''.join(["What country is ", affil, " in? : "]))
    while match_country:
        country = match_country if len(match_country) > 0 else country
        match_country = input(''.join([country,
                        "? (enter for yes or type again): "]))
    if country:
        return(country)
    else:
        return(country_prompt(countries, affil))


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


def save_vocab(which_vocab, voc):
    """
    Function to get the current cntrolled vocabulary for affiliations

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
        json.dump(voc, js)

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
    """
    up_fun = {'affiliation': (parse_affiliations, os.path.join(br_path,
             'biblio_reader', 'map', 'affiliations.json'))}
    geo_dict, new_v = up_fun[which_vocab][0](up_fun[which_vocab][1])
    voc = all_matched_searches(geo_dict, new_v)
    return(voc)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    which_vocab = 'affiliation'
    v = get_vocab(which_vocab)
    v = update_vocab(which_vocab)
    save_vocab(which_vocab, v)
