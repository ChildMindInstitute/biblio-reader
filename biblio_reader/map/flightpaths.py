#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates JSON with affiliation information for mapmaking flightpaths

Authors:
    – Michael Fleischmann, 2017
    – Jon Clucas, 2017

Copyright ©2017, Child Mind Institute (http://childmind.org), Apache v2.0
    License
"""
from itertools import combinations
import json
import os
import sys
br_path = os.path.abspath(os.path.join(__file__,
                          os.pardir, os.pardir, os.pardir))
if br_path not in sys.path:
    sys.path.append(br_path)
from biblio_reader.map.affil_to_json import get_affiliation_json


def main():
    fp = flightpaths(get_affiliation_json()[0])
    with open('flightpaths.json', 'w') as fpfp:
        json.dump(fp, fpfp,
                  sort_keys=True, indent=4, separators=(',', ': ')
                  )


def flightpaths(affiliations):
    """
    Function to find pairs of latitudes and longitude pairs for flightpaths

    Parameter
    ---------
    affiliations: dictionary
        keys are latitude, longitude pairs as strings; values are dictionaries
        of "papers": list of ints (internal paper indices); "affiliations":
        strings as spelled out in papers; "matched searches": terms fed to the
        Google Maps API that returned parent key coordinates

    Returns
    -------
    flightpaths: set
        set of JavaScript Google Maps flightpath strings
    """
    papersets = dict()
    for a in affiliations:
        for p in affiliations[a]['papers']:
            if p in papersets:
                papersets[p].append(a)
            else:
                papersets[p] = [a]
    coords = dict()
    for p in papersets:
        coords[p] = set(combinations(papersets[p], 2))
    flightpaths = dict()
    i = 0
    for p in coords:
        for c in coords[p]:
            s = "".join(["_", str(i)]) if i else ""
            flightpaths["".join(["coauthorship", s])] = [
                                     {'lat': float(c[0].split(',')[0]),
                                      'lng': float(c[0].split(',')[1])},
                                     {'lat': float(c[1].split(',')[0]),
                                      'lng': float(c[1].split(',')[1])}
                                                      ]
            i = i + 1
    return(flightpaths)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
