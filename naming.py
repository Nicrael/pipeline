#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Naming routines
'''

# System modules
import sys
from pathlib import Path
from astropy import log
from astropy.time import Time

ARP = 'arp'
CWD = Path.cwd()


def timestamp(date=False, iso=False):
    '''
    Get a time stamp for strings in iso format.
    If date=True, return only yyyy-mm-dd.
    '''

    time = Time.now()
    string = time.iso[:-4] if not iso else time.isot[:-4]
    if date:
        string = string.split()[0]

    return string


def hist(text=None):
    '''
    Add a history tag of type. Example:
    [ arp.fill_header.newhead ] 2020-04-29T20:06:47 > Created.
    '''
    tsp = timestamp()

    text = 'Updated.' if not text else text
    routine = sys._getframe().f_back.f_code.co_name
    string = f"[ {ARP}.{routine} ] {tsp} > {text}"

    return string


def skeleton(date=False, dry_run=False):
    '''
    Create the directory structure for processed data.
    '''

    tsp = '-'+timestamp(iso=True) if date else ''

    #tsp = timestamp()
    struc = [CWD,
             Path(f'{ARP}-data{tsp}'),
             Path('./reduced'),
             Path('../solved'),
             Path('../master-bias'),
             Path('../master-dark'),
             Path('../master-flat')]

    structure = Path.joinpath(*struc)

    if not dry_run:
        log.info("Creating directory structure.")
        structure.mkdir(parents=True, exist_ok=True)
    else:
        log.warning("Fake creation of directory structure.")
        return

    return structure


def trim(combi):
    '''
    Replacing wired characters from combination of parameters.
    '''
    #tsp = timestamp()

    rep_chars = {'-': '',
                 '"': '',
                 "'": '',
                 ' ': '',
                 '.': '',
                 ':': '=',
                 '{': '',
                 '}': '',
                 '[': '',
                 ']': '',
                 '(': '',
                 ')': '',
                 ',': '.',
                 }

    trimmed_combi = combi
    for k in rep_chars:
        trimmed_combi = str(trimmed_combi).replace(k, rep_chars[k])

    log.debug("Trimming {combi} : {trimmed_combi}.", combi=combi, trimmed_combi=trimmed_combi)
    return trimmed_combi


def output_file(product=None, text=None, counter=None):
    '''
    Creating a default name for output files.
    '''

    if not product:
        product = 'generic'
    product = f'{ARP}.{product}'

    if text:
        text = trim(text)

    if counter:
        text += "."+str(counter).zfill(3)

    ext = '.fits'
    out = f'{product}.{text}{ext}'

    #log.info(f"Creating output file: {output_file}.")
    return out
