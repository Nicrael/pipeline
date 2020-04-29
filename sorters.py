#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
import numpy as np
import itertools as it


from reduction import get_fits_header


def sort(dict_list, keys):
    '''
    Sort a dict list by a list of keys.
    Useful for sorting a table_header by JD
    '''
    sort_by = [keys] if isinstance(keys, str) else keys
    first = lambda x: (tuple(x[s] for s in sort_by) )
    dict_list = sorted(dict_list, key=first)

    return dict_list


def table(pattern, sort_by=None):
    '''
    Take a pattern of fits file names and stacks the headers
    in a dict list. [{FILTER:V, EXPTIME:10}, {FILTER:R, EXPTIME:10}]
    It adds the full path filename.
    Can be sorted by a list of keywords, for example JD.
    '''
    if isinstance(pattern, str):
        pattern = [pattern]


    log.info(f"Stacking {len(pattern)} in a table. It can take a while...")
    heads = [ get_fits_header(f, fast=True) for f in pattern ]
    tabled_heads = [ dict([ [h["name"],h["value"]] for h in H.records()]) for H in heads ]

    log.info("Adding a FULLPATH keyword to retrieve the original file path...")
    full_path = [ {'FULLPATH': f} for f in pattern ]
    tabled_heads = [{**u, **v} for u, v in zip(tabled_heads, full_path)]

    if sort_by:
        tabled_heads = sort(tabled_heads, sort_by)

    return tabled_heads


def group(table_dict, keys, show=None):
    '''
    Take dict list, for example of stacked headers, and extract
    values from a kist of keywords. Can show separated values.
    For example sort unique FILTER keywords and show corresponding JDs.
    '''

    tab = table_dict

    sort_by = [keys] if isinstance(keys, str) else keys
    first = lambda x: (tuple(x[s] for s in sort_by) )
    table_sorted = sorted(tab, key=first)
    result = [ list(g) for k,g in it.groupby(table_sorted, first) ]

    # Avoiding single element tuples
    if isinstance(keys, list) and len(keys) != 1: # ["A","B"]
        values = [ tuple(d.get(k) for k in keys) for d in tab  ]
    elif isinstance(keys, list) and len(keys) == 1: # ["A"]
        values = [ d[keys[0]] for d in tab  ]
    else: # "A"
        values = [ d[keys] for d in tab  ]

    names = show
    if not names:
        res3 = values
    else:
        # Avoiding single element tuples
        if isinstance(names, list) and len(names) != 1: # ["C","D"]
            res2 = [ [ tuple(x[n] for n in names) for x in r] for r in result ]
        elif isinstance(names, list) and len(names) == 1: # ["C"]
            res2 = [ [ x.get(names[0]) for x in r] for r in result ]
        else: # "D"
            res2 = [ [ x.get(names) for x in r] for r in result ]

        res3 = dict(zip(set(values),res2))

    return res3


class minidb():

    def __init__(self, pattern):
        if isinstance(pattern, str):
            pattern = [pattern]
        self.pattern = pattern
        self.heads = [ get_fits_header(f, fast=True) for f in pattern ]
        keys = self.heads[0].keys()
        values = [ [ h.get(k) for h in self.heads ] for k in keys ]
        dic = dict(zip(keys, values))
        dic["ARP FILENAME"] = pattern # adding filename
        self.dic = dic
        self.table = Table(dic) # original
        self.data = self.table
        self.unique = None
        self.names = None

    def group_by(self, keys):
        # if isinstance(keys, str):
        #     keys = [keys]
        self.data = self.table.group_by(keys)
        self.unique = self.data.groups.keys.as_array().tolist()
        return self.data.groups


    def show(self, keys):
        self.names = np.array(self.data[keys]).tolist()
        #self.data = self.table[keys]
        return self.names


class dfits():
    '''
    dfits | fitsort simple clone.
    Uses fast fitsio method by default.
    '''
    def __init__(self, pattern):
        if isinstance(pattern, str):
            pattern = [pattern]
        self.pattern = pattern
        self.heads = [ get_fits_header(f, fast=True) for f in pattern ]
        self.data = self.heads

    def fitsort(self, keys):
        ph = zip(self.pattern,self.heads)
        results = [ (p,(tuple(h[k] for k in keys))) for p,h in ph ]
        self.keys = keys
        self.names = [ r[0] for r in results ]
        self.values = [ r[1] for r in results ]
        self.unique_values = set(self.values)
        self.data = results
        return self

    def unique_names_for(self, value):
        un = [ d[0] for d in self.data if d[1] == value ]
        return un

    def grep(self, value):
        gr = [ d for d in self.data if d[1] == value ]
        return gr
