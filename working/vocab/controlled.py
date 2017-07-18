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
    ll = dict()
    iso_3166_1 = pd.read_csv(os.path.abspath(os.path.join(__file__, os.pardir,
                 "ISO_3166_1.csv")), na_filter=False)
    iso_3166_2_us = pd.read_csv(os.path.abspath(os.path.join(__file__,
                    os.pardir, "ISO_3166_2_US.csv")), na_filter=False)
    doublecheck = {"India", "Georgia", "Mexico", "SRB", "BIH", "CA", "Jersey",
                   "Israel"}\
                  | {country['Alpha-2 code'] for country in iso_3166_1.to_dict(
                  orient='records')} | {country['Alpha-3 code'] for country in
                  iso_3166_1.to_dict(orient='records')}
    singlecheck = {"".join([", ", country]) for country in doublecheck if (
                  country not in {country['Alpha-2 code'] for country in
                  iso_3166_1.to_dict(orient='records')} and country not in
                  {country['Alpha-3 code'] for country in iso_3166_1.to_dict(
                  orient='records')})} | {"China", "Taiwan", "PRC"} | {state[
                  'Code'] for state in iso_3166_2_us.to_dict(orient='records')}
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
                '?': 'unknown', 'Taiwan': 'TW'}
    us = {'US', 'USA', 'United States', 'U.S.A'}
    us_states = {state['Code'][-2:]: state['Code'] for state in
                iso_3166_2_us.to_dict(orient='records')}
    for state in us_states:
        if state not in countries:
            countries[state] = us_states[state]
    country_count = {**{country['Alpha-2 code']: 0 for country in
                    iso_3166_1.to_dict(orient='records')}, **{state['Code']: 0
                    for state in iso_3166_2_us.to_dict(orient='records')}}
    for country in countries:
        for keycountry in countries:
            if country in keycountry and not country == keycountry and country\
               not in singlecheck:
                if len(country) > 3 and country not in doublecheck:
                   singlecheck.add(country)
                else:
                    doublecheck.add(country)
                    singlecheck.add("".join([", ", country]))
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
                    latlong[affil] = {"country": country_prompt(countries,
                                     iso_dict, affil, doublecheck, singlecheck,
                                     us_states)}
    for i, affil in enumerate(latlong):
        matched_country = None
        canon_match = None
        n = []
        if not matched_country and affil and latlong[affil] and "country" not  \
           in latlong[affil]:
            while not matched_country:
                matched_country = country_prompt(countries, iso_dict, affil,
                                  doublecheck, singlecheck, us_states)
                for keycountry in countries:
                    if matched_country in keycountry and not matched_country  \
                       == keycountry:
                        doublecheck.add(matched_country)
                        singlecheck.add("".join([", ", matched_country]))
                n.append(matched_country)
            if len(n) > 1:
                matched_country = None
                for m in n:
                    matched_country = country_prompt(countries, iso_dict,
                                      m, doublecheck, singlecheck, us_states) \
                                      if not matched_country else             \
                                      matched_country
            else:
                matched_country = n[0]
            if matched_country in countries and affil in latlong:
                canon_match = countries[matched_country]
                latlong[affil]["country"] = canon_match
            elif matched_country in countries:
                latlong[affil] = dict()
                latlong[affil]["country"] = canon_match
            else:
                if(len(countries)) and affil not in latlong:
                    print("\n")
                    for cc in sorted(list({cc for c, cc in countries.items(
                                     )})):
                        print(iso_dict[cc][0], end="; ")
                    print("\n")
                    print("If that country is the same as one of the "
                          "above, please enter the match. Otherwise, just "
                          "press `enter`")
                    canon_match = countries[input(
                                  "(please type the full name of the "
                                  "match, matching case): ")]
                canon_match = matched_country if not canon_match else     \
                              canon_match
                if matched_country not in affil:
                    matched_country = input("How's that spelled here? ")
                if not matched_country or len(matched_country) == 0:
                    matched_country = canon_match
                latlong[affil]["country"] = canon_match
            if canon_match:
                countries[matched_country] = canon_match
                latlong[affil]["country"] = canon_match
                if canon_match not in countries:
                    countries[canon_match] = canon_match
            if canon_match in country_count:
                country_count[canon_match] += 1
            else:
                print(matched_country)
        elif "country" in latlong[affil]:
            canon_match = latlong[affil]["country"]
            if canon_match in country_count:
                country_count[canon_match] += 1
                if "country" not in latlong[affil]:
                    latlong[affil]["country"] = canon_match
        if canon_match in country_count:
            print(": ".join([iso_dict[canon_match][0], str(country_count[
                     canon_match]), affil]))
        else:
            print(affil)
        if "country" not in latlong[affil]:
            latlong[affil]["country"] = canon_match

    for coord, v in affiliations.items():
        for term in v["matched searches"]:
            for match, lv in latlong.items():
                if "country" not in affiliations[coord] and term in lv['match']:
                    affiliations[coord]["country"] = lv["country"]
    save_heatmap_data(affiliations)
    return(affiliations, country_count)


def country_prompt(countries, iso_dict, affil, doublecheck=["India", "CA",
                   "MA", "Georgia", "Mexico", "IN", "PR"], singlecheck=[],
                   states={}, us={'US', 'USA', 'United States', 'U.S.A'}):
    """
    Function to prompt user if affiliation is in a given country, and, if not,
    to identify the affiliation.

    Parameters
    ----------
    countries: dictionary
        countries vocabulary

    iso_dict: dictionary
        spelled-out country names

    affil: string
        affiliation

    doublecheck: list of strings
        countries to doublecheck if found

    doublecheck: list of strings
        strings to not doublecheck if found

    states: dictionary
        dictionary of state abbreviatons

    Returns
    -------
    country: string or None
        confirmed country
    """
    for country in countries:
        if country and "".join([" ", country]) in affil:
            for u in us:
                if u in affil and country in affil:
                    if affil.endswith(u) and affil not in ["VA", "PR"]:
                        return("".join(["US-", countries[country][-2:]]))
                    return(state_prompt(countries, iso_dict, affil, states))
            if country not in doublecheck or (country in singlecheck and
               "Indiana" not in affil):
                if countries[country] in us:
                    return(state_prompt(countries, iso_dict, affil, states))
                elif "Indiana" in affil:
                    return("US-IN")
                else:
                    return(country)
            if affil.endswith(country) and len(country) > 2:
                return(country)
            match_country = "maybe"
            while(len(match_country) > 0 and match_country.upper()
                  not in ["Y", "N"]):
                if country == affil[-len(country):] and country in countries:
                    return(countries[country])
                match_country = input(''.join(["Is ", affil, " in the state",
                                " \"", iso_dict[countries[country]][0],
                                "\"? (Y / N) : "]))
            if len(match_country) == 0 or match_country.upper() == "Y":
                if countries[country] in us:
                    return(state_prompt(countries, iso_dict, affil, states))
                return(country)
            if countries[country] not in us and ((country and country[0:2] ==
               "US" and "".join([" ", country[-2:]]) in affil) or ("".join([
               " ", country]) in affil)) and "".join(["US-", countries[country
               ][-2:]]) in countries:
               for u in us:
                   if ", ".join([countries[country][-2:], u]) in affil:
                       return("".join(["US-", countries[country][-2:]]))
               if "US" in iso_dict[countries[country]] or input(''.join([
                  "Is ", affil, " in the state \"", iso_dict["US"][0],
                  "\"? (Y / N) : "])).upper() != 'N':
                   return(state_prompt(countries, iso_dict, affil, states))
    match_country = None
    while not(match_country):
        match_country = input(''.join(["What state is ", affil, " in? : "]))
    while match_country:
        country = match_country if len(match_country) > 0 else country
        match_country = input(''.join([country,
                        "? (enter for yes or type again): "]))
    if country:
        return(country)
    else:
        return(country_prompt(countries, iso_dict, affil, doublecheck,
               singlecheck, states))


def state_prompt(countries, iso_dict, affil, states, us={'US', 'USA',
                 'United States', 'U.S.A'}):
    """
    Function to prompt user if affiliation is in a given US state, and, if not,
    to identify the affiliation.

    Parameters
    ----------
    countries: dictionary
        countries vocabulary

    iso_dict: dictionary
        spelled-out state names

    affil: string
        affiliation

    states: dictionary
        dictionary of state abbreviatons

    Returns
    -------
    country: string or None
        confirmed state
    """
    for country in countries:
        for u in us:
            if u in affil and country in affil:
                if "".join(["US-", countries[country][-2:]]) in countries and  \
                   country not in ["VA", "PR"]:
                       return("".join(["US-", countries[country][-2:]]))
        if country not in us and countries[country] != "US" and ((country and
           country[0:2] == "US" and "".join([" ", country[-2:]]) in affil) or (
           "".join([" ", country]) in affil)) and "".join(["US-", countries[
           country][-2:]]) in countries:
            match_country = "maybe"
            while(len(match_country) > 0 and match_country.upper() not in ["Y",
                  "N"]):
                country = "".join(["US-", country]) if "".join(["US-", country]
                          ) in iso_dict else country
                if country == affil[-len(country):] and country in countries:
                    return(countries[country])
                match_country = input(''.join(["Is ", affil,
                                " in the US state \"", iso_dict[countries[
                                country]][0], "\"? (Y / N) : "]))
                if len(match_country) == 0 or match_country.upper() == "Y":
                    return("".join(["US-", countries[country][-2:]]))
    match_country = None
    while not(match_country):
        match_country = input(''.join(["What US state is ", affil, " in? : "]))
    while match_country:
        country = match_country if len(match_country) > 0 else country
        match_country = input(''.join([country,
                        "? (enter for yes or type again): "]))
    if country and country in countries and "".join(["US-", countries[country][
       -2:]]) in countries:
        return("".join(["US-", countries[country][-2:]]))
    elif country:
        return(country)
    else:
        return(None)


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
