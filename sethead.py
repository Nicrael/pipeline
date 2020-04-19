#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
import json

# Our modules
from reduction import get_fits_header, is_number, to_number

def sethead(head, key, val):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file
    
    In [23]: "Hello %.2f" % 123.123
    Out[23]: 'Hello 123.12'
    
    In [30]: 'Hello, {:.2f}'.format(123.123)
    Out[30]: 'Hello, 123.12'
    '''

    with open('cerbero-merged-array.json') as jfile:
        head_format = json.load(jfile)['primary']
    
    card = [c for c in head_format if c["name"] == key]
        
    if not card:
        form = '{}'
        comm = None
    else:
        form = '{'+card[0]["format"]+'}'
        comm = card[0]["comment"]
    
    if is_number(val):
        value = form.format(val)
        value = float(value)
    else:
        value = val

    #print(key, val, "→", form, "→", value)
      
    if not card:
        head[key] = value
    else:
        head[key] = (value, card[0]["comment"])

    return head
