#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sort_by.py

A function to sort objects.json and save the sorted object file for easier
manual editing of said file.

Author:
    –Jon Clucas, 2017
    
   © 2016, Child Mind Institute, The MIT License
"""
import json, operator

def sort_oj():
    """
    A function to sort objects.json and save the sorted file.
    """
    sort_keys = {"method" : 1,
                 "software" : 2,
                 "paper" : 3,
                 "protocol": 4,
                 "database": 5,
                 "domain" : 6}
    with open("objects.json", "r") as oj:
        data = json.load(oj)
    new_json = sorted(data, key=operator.itemgetter('name'))
    data = sorted(new_json, key=lambda x: sort_keys[x['type']])
    with open('objects.json', 'w') as f:
        json.dump(data, f, sort_keys=True, indent=4, separators=(',', ':'))

# ============================================================================
if __name__ == '__main__':
    sort_oj()