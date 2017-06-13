#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 17:09:49 2017

@author: jon.clucas
"""

import json

with open('lpt.json', 'r') as f:
    data = json.load(f)

new_json = []

for node in data["_nodes"]:
    new_json.append({"name" : node, "type": "unclassified", "depends" : []})

for suc in data["_sucs"]:
    for key in data["_sucs"][suc].keys():
        for item in new_json:
            if (item['name'] == key):
                item['depends'].append(suc)

with open('new_objects.json', 'w') as f:
    json.dump(new_json, f)