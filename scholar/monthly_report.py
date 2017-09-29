#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
monthly_report.py

Created on Thu Sep 21 09:28:06 2017

@author: jon.clucas
"""
from datetime import date
import subprocess

global keywords
keywords = [
               "fcon_1000.projects.nitrc.org",
               "Rockland Sample",
               "1000 Functional Connectomes",
               "International Neuroimaging Data-Sharing Initiative",
               "Autism Brain Imaging Data Exchange",
               "ADHD-200",
               "Consortium for Reproducibility and Reliability"
           ]

command = "python biblio_reader.py -s \"{0}\" --before {1} --after {2} -o \""\
          "{3}\"".format(
              ", ".join([k for k in keywords]),
              str(date.today().year),
              str(date.today().year),
              str(date.today())
          )
print(command)
subprocess.call(command, shell=True)
