#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy import log
from astropy.time import Time
from pathlib import Path
import sys

arp = 'arp'
cwd = Path.cwd()

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
    ts = timestamp()

    text = 'Updated.' if not text else text
    routine = sys._getframe().f_back.f_code.co_name
    string = f"[ {arp}.{routine} ] {ts} > {text}"

    return string


def skeleton(date=False, dry_run=False):
    '''
    Create the directory structure for processed data.
    '''

    ts = '-'+timestamp(iso=True) if date else ''

    #ts = timestamp()
    struc = [ cwd,
              Path(f'{arp}-data{ts}'),
              Path('./reduced'),
              Path('../solved'),
              Path('../master-bias'),
              Path('../master-dark'),
              Path('../master-flat') ]

    structure = Path.joinpath(*struc)

    if not dry_run:
        log.info("Creating directory structure.")
        structure.mkdir(parents=True,exist_ok=True)
    else:
        log.warn("Fake creation of directory structure.")
        return

    return structure


def trim(combi):
    '''
    Replacing wired characters from combination of parameters.
    '''
    ts = timestamp()

    rep_chars = {'-' : '',
                 '"' : '',
                 "'" : '',
                 ' ' : '',
                 '.' : '',
                 ':' : '=',
                 '{' : '',
                 '}' : '',
                 '[' : '',
                 ']' : '',
                 '(' : '',
                 ')' : '',
                 ',' : '.',
                 }

    trimmed_combi = combi
    for k in rep_chars.keys():
        trimmed_combi = str(trimmed_combi).replace(k,rep_chars[k])

    log.debug(f"Trimming {combi} : {trimmed_combi}.")
    return trimmed_combi


def output_file(product=None, text=None, counter=None):
    '''
    Creating a default name for output files.
    '''

    if not product:
        product = 'generic'
    product = f'{arp}.{product}'

    if text:
        text= trim(text)

    if counter:
        text += "."+str(counter).zfill(3)

    ext = '.fits'
    output_file = f'{product}.{text}{ext}'

    #log.info(f"Creating output file: {output_file}.")
    return output_file
