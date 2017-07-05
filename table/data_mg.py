#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Used to isolate specific items from the Pandas data structure
Can be isolated based on publication type, whether it uses the dataset, the
    affiliation, author, etc.


Authors:
    – Michael Fleischmann, 2017
    – Jon Clucas, 2017

Copyright 2017, Child Mind Institute (http://childmind.org), Apache v2.0
    License
"""

import os
import sys
br_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
if br_path not in sys.path:
    sys.path.append(br_path)


def table_data(data, data_use, journal_category, cols, sort_cols, write=True,
               fn='Table_data.csv'):
    """
    Parameters
    ----------
    data: pandas dataframe
        master data

    data_use: string
        ∈ ['Y', 'S', 'N', 'I']

    journal_category: string
        "Journal", "Thesis", etc.

    cols: list of strings
        column headers to keep in subset

    sort_cols: list of strings
        column headers for columns to sort by

    write: boolean
        save csv? default=True

    fn: string
        output filename, default="Table.csv"
    """
    data = data[(data['Data Use'] == data_use) & (data['Journal Category'] ==
           journal_category)]
    table_data = data[cols].sort_values(sort_cols)
    if write:
        table_data.to_csv(path_or_buf=fn, index=False)
    return(table_data)


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import manager as mg
    table_data(mg.get_data(), "Y", "Thesis", ['Title', 'Authors', 'Year',
               'URL'], ['Year', 'Title'])