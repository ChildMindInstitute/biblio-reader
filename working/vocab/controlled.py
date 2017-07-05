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
br_path = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(
          __file__))))
if br_path not in sys.path:
    sys.path.append(br_path)
import json
import manager as mg
import biblio_reader.map.affil_to_json as atj


def all_matched_searches(affiliations, de_facto_affiliations):
    """
    Function to create a dictionary for all matched geolocation search terms.

    Parameters
    ----------
    affiliations: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates

    de_facto_affiliations: set
        set of strings of affiliations

    Returns
    -------
    latlong: dictionary
        dictionary of affiliation, latitude-longitude string pairs
    """
    ll = dict()
    for coord, v in affiliations.items():
        for term in v["matched searches"]:
            ll[term] = coord
    latlong = dict()
    for affil in de_facto_affiliations:
        if affil in ll:
            latlong[affil] = ll[affil]
        else:
            aff_l, aff_r = parse_affiliation(affil)
            if aff_r in ll:
                latlong[affil] = {aff_r: ll[aff_r]}
            else:
                while(aff_r not in ll and len(aff_r.split(", ")) > 1):
                    aff_l, aff_r = parse_affiliation(aff_r)
                if(aff_r in ll):
                    latlong[affil] = {aff_r: ll[aff_r]}
                else:
                    latlong[affil] = None
    return(latlong)


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
    v_path = os.path.dirname(path)
    if not os.path.exists(v_path):
        os.makedirs(v_path)
    with open(path, 'w') as js:
        json.dump(voc, js)


def update_vocab(which_vocab, voc):
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
    voc: dictionary
        dictionary with de facto spelling string keys and canonical spelling
        string values, updated
    """
    up_fun = {'affiliation': (parse_affiliations, os.path.join(br_path,
             'biblio_reader', 'map', 'affiliations.json'))}
    geo_dict, new_v = up_fun[which_vocab][0](up_fun[which_vocab][1])
    matches = all_matched_searches(geo_dict, new_v)
    for term in new_v:
        if term not in voc:
            if term in matches:
                voc[term] = matches[term]
            else:
                voc[term] = term
    return(voc)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    which_vocab = 'affiliation'
    v = get_vocab(which_vocab)
    v = update_vocab(which_vocab, v)
    save_vocab(which_vocab, v)
