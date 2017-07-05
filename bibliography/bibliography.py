#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script creates a human-readable html bibliography from a master data csv

Author:
    – Jon Clucas, 2017

Copyright ©2017, Child Mind Institute (http://childmind.org), Apache v2.0
    License
"""

import argparse
import os
import pandas as pd
import sys
br_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
if br_path not in sys.path:
    sys.path.append(br_path)
from table.data_mg import table_data
from manager import get_data

def apa_format(auth, year, small, large, url="", highlight=None):
    """
    Function to format data into APA-ish format.

    Parameters
    ----------
    auth: string
        authors

    year: string
        year

    small: string
        small work

    large: string
        large work

    url: string
        web location

    highlight: string or None
        HTML color code for highlight color

    Returns
    -------
    apaish: string
        formatted html string
    """
    hightlight = "" if not highlight else "".join([' style="background-color:',
                 highlight, '"'])
    url = "".join([' [<a href="', url, '">link</a>].']) if len(url) > 0 else ""
    apaish = "".join(['''<p class="apaish"''', hightlight, """>""",
             auth.rstrip("."), ". (20", year[-2:] + "). "])
    if large and len(large) > 0:
        apaish = "".join([apaish, '"', small.rstrip("."), '." <cite>',
                 large.rstrip("."), "</cite>.", url, "</p>"])
    else:
        apaish = "".join([apaish, '<cite>', small.rstrip("."), "</cite>.", url,
                 "</p>"])
    return(apaish)


def build_html(data, html_out, cats=None, uses=['Y', 'S', 'N'], sets=None,
               title="Bibliography"):
    """
    Function to create an HTML page.

    Parameters
    ----------
    data: pandas dataframe
        dataframe to build from

    html_out: string
        path to save HTML to

    cats: list of strings or set of strings or None
        journal categories to include. If none, include all in no particular
        order

    sets: list of strings
        datasets

    uses: list of strings
        data uses to include from ['Y', 'S', 'N', 'I'], default=['Y', 'S', 'N']

    title: string
        title of document, default="Bibliography"

    Returns
    -------
    html_open: string
        opening html
    """
    cats = {cat for cat in data['Journal Category']} if not cats else cats
    html_string = html_open(title)
    for cat in cats:
        for use in uses:
            html_string = html_string + html_bib_block(data, use, cat, sets)
    html_string = html_string + html_close()
    with open(html_out, 'w') as h:
        h.write(html_string)


def html_bib_block(data, data_use, journal_category, sets=None):
    """
    Function to create a block of HTML bibliography.

    Parameters
    ----------
    data: pandas dataframe
        dataframe to build from

    data_use: string
        ∈ ['Y', 'S', 'N', 'I']

    journal_category: string
        "Journal", "Thesis", etc.

    sets: list of strings or None
        data sets, optional

    Returns
    -------
    html_block: string
        html for bibliography section
    """
    html_string = ""
    if sets and len(sets) > 1:
        for set_name in sets:
            html_string += (html_bib_block(data, data_use, journal_category, [
                           set_name]))
    category = {"Journal": "Peer-reviewed journal articles",
                "Proceeding": "Articles in conference proceedings",
                "Preprint": "Preprint articles",
                "Thesis": "Theses and dissertations"}
    set_names = {'FCP': '1,000 Functional Connectomes',
                 'ADHD200': 'ADHD-200',
                 'NKI': 'NKI-Rockland',
                 '': 'unspecified'}
    set_name = 'our' if not sets else set_names[sets[0]] if sets[0] in        \
                set_names else sets[0]
    suffix = {'Y': ' '.join(['that used', set_name, 'data']),
              'N': ' '.join(['that cited', set_name, 'data']),
              'S': ' '.join(['that used', set_name, 'scripts'])}
    bib_data = table_data(data, data_use, journal_category, ['Authors', 'Year',
               'Title', 'Journal', 'Sets', 'URL', 'Contributor'], ['Authors',
               'Year', 'Title', 'Journal'],
               False).fillna("")
    if bib_data.shape[0] == 0:
        return("")
    journal_category = category[journal_category] if journal_category in      \
                       category else journal_category
    data_use = suffix[data_use] if data_use in suffix else ""
    html_string += "".join(["""
                <p><strong>""", journal_category, " ", data_use,
                """</strong>"""])
    for i, row in bib_data.iterrows():
        if not sets or (len(sets[0]) and sets[0] in row.loc['Sets']) or sets[0
                       ] == row.loc['Sets']:
            q = "#f2cd32" if ('others' in row.loc['Authors'] or set_name ==
                "unspecified" or "books.google.com" in row.loc['URL']) else   \
                "#00c1d5" if row.loc['Contributor'] == 'Contributor' else None
            html_string += apa_format(str(row.loc['Authors']), str(row.loc[
                           'Year']), str(row.loc['Title']), str(row.loc[
                           'Journal']), row.loc['URL'], q)
    html_string += """
                </p>"""
    return(html_string)


def html_close():
    """
    Function to create the tail of an HTML page.

    Parameters
    ----------
    None

    Returns
    -------
    html_close: string
        close html
    """
    return("""

        </section>

    </div>

</main>
</body>
</html>""")


def html_open(title="Bibliography", description=""):
    """
    Function to create the head of an HTML page.

    Parameters
    ----------
    title: string
        page title, default="Bibliography"

    description: string
        page description, default=""

    Returns
    -------
    html_open: string
        opening html
    """
    description_p = "".join(["""<p class="col-head__intro">""", description,  \
                    """</p><!-- /.intro -->"""]) if len(description) > 0 else \
                    description
    return("".join(["""<!DOCTYPE html>
<!--[if lte IE 8 ]>
<html lang="en" class="no-js oldie">
<![endif]-->
<!--[if IE 9 ]>
<html lang="en" class="no-js ie9">
<![endif]-->
<!--[if (gt IE 9)|!(IE)]><!-->
<html lang="en" class="no-js">
<!--<![endif]-->
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
    <meta http-equiv="X-UA-Compatible" content="IE=edge,chrome=1">

    <link rel="apple-touch-icon-precomposed" sizes="152x152" href="https://"""
"""27c2s3mdcxk2qzutg1z8oa91-wpengine.netdna-ssl.com/wp-content/themes/"""
"""childmind/assets/img/favicons/apple-touch-icon-152x152.png" />
    <meta name="msapplication-TileColor" content="#FFFFFF" />
    <meta name="msapplication-TileImage" content="https://childmind.org/"""
"""wp-content/themes/childmind/assets/img/favicons/mstile-144x144.png" />
<title>""", title, '''</title>
<meta name="description" content="''', description, '''"/>
<meta property="og:locale" content="en_US" />
<meta property="og:type" content="website" />
<meta property="og:title" content="Facebook Logo" />
<meta property="og:description" content="Child Mind Institute" />
<meta property="og:url" content="https://childmind.org" />
<meta property="og:site_name" content="Child Mind Institute" />
<meta property="og:image" content="https://childmind.org/wp-content/uploads/'''
'''cma_vote_social_img_share.png" />
<meta name="twitter:card" content="summary" />
<meta name="twitter:description" content="''', description, '''" />
<meta name="twitter:title" content="''', title, '''" />
<meta name="twitter:site" content="@ChildMindDotOrg" />
<meta name="twitter:image" content="https://childmind.org/wp-content/'''
'''uploads/cma_vote_social_img_share.png" />

<link rel='stylesheet' id='clicktofbshare-css'  href='https://'''
"""27c2s3mdcxk2qzutg1z8oa91-wpengine.netdna-ssl.com/wp-content/plugins/click"""
"""-to-fbshare/assets/css/styles.css?ver=4.4.8' type='text/css' media='all' />
<link rel='stylesheet' id='tm_clicktotweet-css'  href='https://"""
"""27c2s3mdcxk2qzutg1z8oa91-wpengine.netdna-ssl.com/wp-content/plugins/click"""
"""-to-tweet-by-todaymade/assets/css/styles.css?ver=4.4.8' type='text/css' """
"""media='all' />
<link rel='stylesheet' id='dashicons-css'  href='https://"""
"""27c2s3mdcxk2qzutg1z8oa91-wpengine.netdna-ssl.com/wp-includes/css/"""
"""dashicons.min.css?ver=4.4.8' type='text/css' media='all' />
<link rel='stylesheet' id='site-css-css'  href='https://"""
"""27c2s3mdcxk2qzutg1z8oa91-wpengine.netdna-ssl.com/wp-content/themes/"""
"""childmind/style.css?ver=4.4.8' type='text/css' media='all' />
<script type='text/javascript' src='https://27c2s3mdcxk2qzutg1z8oa91-"""
"""wpengine.netdna-ssl.com/wp-content/themes/childmind/src/js/vendor/"""
"""jquery.js?ver=1.11.3'></script>
<script type='text/javascript' src='https://27c2s3mdcxk2qzutg1z8oa91-"""
"""wpengine.netdna-ssl.com/wp-content/themes/childmind/assets/js/modernizr"""
""".js?ver=3.0.0'></script>
<link rel='https://api.w.org/' href='https://childmind.org/wp-json/' />
<style>
a {
    color: #242a6a;
}

body {
    background-color: #efefef;
}

cite {
    font-style: italic;
}

strong {
    font-size:20px;
    font-size:1.1875rem;
    font-weight: bold;
    line-height: 2.5;
}

.apaish {
    padding-left: 4em;
    text-indent: -4em;
}
</style>
</head>

<body id="body" class="home blog facetwp-is-loading">
<main class="layout--b">

    <div class="row">

        <section>

            <div>

                <h1 class="title">""", title, """</h1><!-- /.title -->""",
                description_p, """
           </div>
"""]))


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='data_in & html_out')
    parser.add_argument('data_in', type=str, nargs='?', default=None,
                        help="path to master CSV (optional)")
    parser.add_argument('html_out', type=str, nargs='?', default=
                        os.path.join(br_path, 'bibliography', 'index.html'),
                        help="save path for HTML file (optional)")
    arg = parser.parse_args()
    if arg.data_in:
        data = pd.read_csv(arg.data_in)
    else:
        data = get_data()
    build_html(data, arg.html_out, ["Journal", "Proceeding", "Thesis",
               "Preprint"], sets=["FCP", "ABIDE", "ADHD200", "NKI", "CORR", ""],
               title="Forward citations")
    build_html(data, os.path.join(os.path.dirname(arg.html_out), "extras.html"
               ), uses=["I"], title="Irrelevant & duplicate references")
