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


def apa_format(auth, year, small, large, url=""):
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

    Returns
    -------
    apaish: string
        formatted html string
    """
    apaish = "".join(["""<p class="apaish">""", auth.rstrip("."), ". (20",
             year[-2:] + "). "])
    if large and len(large) > 0:
        apaish = "".join([apaish, '"', small.rstrip("."), '." <cite>',
                 large.rstrip("."), "</cite>.</p>"])
    else:
        apaish = "".join([apaish, '<cite>', small.rstrip("."), "</cite>.</p>"])
    return(apaish)


def build_html(data, html_out):
    """
    Function to create an HTML page.

    Parameters
    ----------
    data: pandas dataframe
        dataframe to build from

    html_out: string
        path to save HTML to

    Returns
    -------
    html_open: string
        opening html
    """
    html_string = (html_open() +
                  html_bib_block(data, 'Y', 'Journal') +
                  html_close())
    with open(html_out, 'w') as h:
        h.write(html_string)


def html_bib_block(data, data_use, journal_category):
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

    Returns
    -------
    html_block: string
        html for bibliography section
    """
    category = {"Journal": "Peer-reviewed journal articles"}
    suffix = {'Y': 'that used our data'}

    bib_data = table_data(data, data_use, journal_category, ['Authors', 'Year',
               'Title', 'Journal'], ['Authors', 'Year', 'Title', 'Journal'],
               False).fillna("")

    html_string = "".join(["""
                <p><strong> """, category[journal_category], " ", suffix[
                  data_use], """</strong>"""])
    for i, row in bib_data.iterrows():
        html_string += apa_format(str(row.loc['Authors']), str(row.loc['Year']
                       ), str(row.loc['Title']), str(row.loc['Journal']))
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
cite {
    font-style: italic;
}

strong {
    font-weight: bold;
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

        <section class="main-content" itemprop="articleBody">

            <div class="col-head">

                <h1 class="title">""", title, """</h1><!-- /.title -->

                <p class="col-head__intro">""",
           description, """</p><!-- /.intro -->

           </div>
"""]))


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='data_in & html_out')
    parser.add_argument('data_in', type=str, nargs='?', default=None, help=
                        "path to master CSV (optional)")
    parser.add_argument('html_out', type=str, nargs='?', default=os.path.join(
                        br_path, 'bibliography', 'index.html'), help=
                        "save path for HTML file (optional)")
    arg = parser.parse_args()
    if arg.data_in:
        data = pd.read_csv(arg.data_in)
    else:
        data = get_data()
    build_html(data, arg.html_out)
