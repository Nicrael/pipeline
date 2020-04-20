#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy.io import fits
import json


def sethead(head, key, val):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file

    In [30]: 'Hello, {:.2f}'.format(123.129)
    Out[30]: 'Hello, 123.13'

    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-array.json') as jfile:
        head_format = json.load(jfile)['primary']

    card = [c for c in head_format if c["name"] == key]
    form = '' if not card else card[0]["format"]

    if 'd' in form:
        value = int(round(float(val)))
    elif 'f' in form:
        value = round(val, [ int(f) for f in form if f.isdigit() ][0])
    else: # 's' in form:
        value = f'{val:{form}}'

    head[key] = value if not card else (value, card[0]["comment"])

    print(key, val, "→", form, "→", value)


    return head
