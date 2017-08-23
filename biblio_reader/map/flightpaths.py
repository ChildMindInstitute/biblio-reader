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
import os
import sys
br_path = os.path.abspath(os.path.join(__file__,
                          os.pardir, os.pardir, os.pardir))
if br_path not in sys.path:
    sys.path.append(br_path)
from biblio_reader.map.affil_to_json import get_affiliation_json


def main():
    fp = flightpaths(get_affiliation_json()[0])
    with open("map_source.js", 'r') as map_source_js:
        map_source = map_source_js.read()
        if "// Flightpaths begin" in map_source:
            map_source = "".join([map_source.split("// Flightpaths begin")[0],
                                  "// Flightpaths begin\n",
                                  "".join(fp),
                                  "// Flightpaths end\n",
                                  map_source.split("// Flightpaths end\n")[1]])
        else:
            print("\n".join(["Please mark the flightpaths section with ",
                             "// Flightpaths begin",
                             "and",
                             "// Flightpaths end",
                             "comments in map_source.js"]))
    with open("map_source.js", 'w') as map_source_js:
        map_source_js.write(map_source)


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
        set of kavascript Google Maps flightpath strings
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
    flightpaths = set()
    i = 0
    for p in coords:
        for c in coords[p]:
            s = str(i) if i else ""
            flightpaths.add("".join([
                            "  var flightPlanCoordinates", s, " = [\n",
                            "    {lat: ", c[0].split(',')[0],
                            ", lng: ", c[0].split(',')[1], "},\n",
                            "    {lat: ", c[1].split(',')[0],
                            ", lng: ", c[1].split(',')[1], "}\n",
                            "  ];\n\n",
                            "  var flightPath", s,
                            " = new google.maps.Polyline({\n",
                            "    path: flightPlanCoordinates", s, ",\n",
                            "    geodesic: true,\n",
                            "    strokeColor: '#f9e28a',\n",
                            "    strokeOpacity: 0.5,\n",
                            "    strokeWeight: 2\n",
                            "    });\n\n",
                            "    flightPath", s, ".setMap(map);\n\n"]))
            i = i + 1
    return(flightpaths)

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
