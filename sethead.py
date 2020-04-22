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
        head_format = json.load(jfile)# ['primary']

    card = [fits.Card(**c) for c in head_format if key in c['keyword']][0]
    #form = '' if not card else card[0]["format"]

    print(key)
    if card:
        if 'd' in card.value:
            card.value = int(round(float(val)))
        elif 'f' in card.value:
            card.value = round(val, [ int(v) for v in card.value if v.isdigit() ][0])
        else: # 's' in form:
            card.value = f'{val:{card.value}}'
    else:
         print(f'{key} not in dict')
     #head[key] = value if not card else (value, card[0]["comment"])

    #print(key, val, "→", form, "→", value)


    return head
