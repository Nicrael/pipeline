#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules

from astropy.io import fits, ascii
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import ccdproc as ccdp
import fitsio
import numpy as np

##########################################################################

class dfits():
    '''
    dfits | fitsort simple clone.
    Uses fast fitsio method by default.
    '''
    def __init__(self, pattern):
        self.pattern = pattern
        self.heads = [ get_fits_header(f, fast=True) for f in pattern ]
        #self.heads = [ fitsio.read_header(f, choose_hdu2(f)) for f in pattern ]

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


##########################################################################


def choose_hdu(filename, fast=False):
    '''
    Detect whether the fits file is compressed with
    fpack, and choose the right HDU.
    fast: Alternative mode based on fitsio
    '''
    if fast:
        finfo = fitsio.FITS(filename) # Object
        finfo_list = [ f.get_extnum() for f in finfo if f.is_compressed() ]
    else:
        finfo = fits.info(filename, output=False) # List of tuples.
        finfo_list = [ f[0] for f in finfo if 'COMPRESSED_IMAGE' in f ]

    if not finfo_list:
        return 0 # finfo # 0 if not compressed
    else:
        return 1 # finfo_list[0] # 1 if compressed


def get_fits_header(filename, fast=False):
    '''
    Return the header of the fits file.
    fast uses fitsio.
    '''
    which_hdu = choose_hdu(filename, fast=fast)
    if fast:
        #header = fitsio.read_header(filename, which_hdu)
        with fitsio.FITS(filename) as f:
            header = f[which_hdu].read_header()
    else:
        header = fits.getheader(filename, which_hdu)
    #print(f"Getting header from {filename}")
    return header


def get_fits_data(filename, fast=False):
    '''
    Return the data of the fits file.
    If fitsio=True, use fitsio.
    '''
    which_hdu = choose_hdu(filename, fast=fast)
    if fast:
        #data = fitsio.read(filename, which_hdu)
        with fitsio.FITS(filename) as f:
            data = f[which_hdu].read()
    else:
        data = fits.getdata(filename, which_hdu)
    #print(f"Getting data from {filename}")
    return data


def is_keyval_in_header(head, key, val):
    '''
    Alternative for 'x' in header and header['x'] == 'y'
    '''
    return key in head and head[key] == val


def is_keyval_in_file(filename, key, val):
    '''
    Alternative for 'x' in fits and fits['x'] == 'y'
    '''
    header = get_fits_header(filename)
    return is_keyval_in_header(header, keyword, value)


def write_fits(data, filename, header=None):
    '''
    Write a fits file.
    It adds a checksum keyword.
    '''

    if header is None:
        hdu = fits.PrimaryHDU(data)
    else:
        hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(filename, overwrite=True, checksum=True)
    return hdu


def mask(data, sigma=3, output_file=None, header=None):
    '''
    Create a bad pixel mask
    '''
    mask = sigma_clip(data, masked=True).mask.astype(int)
    if output_file:
        write_fits(mask, output_file, header=header)

    return mask


def mask_reg(data, sigma=3, output_file=None):
    '''
    Create a bad pixel region table
    '''

    y,x = np.where(data == True)
    p = np.repeat("point ", y.size)
    t = [p, x+1,y+1]

    table = Table(t, names=['# ','## ','###']) # bleah

    if output_file:
        ascii.write(table, output_file, overwrite=True)

    return table


def combine(keys, pattern, prod, method='average', normalize=False, fast=True):
    '''
    Combine a pattern of images using average (default) or median.
    Loops over a list of keywords. normalize=True to combine flats.
    '''

    print('Getting headers of all files in pattern')
    dlist = dfits(pattern).fitsort(keys)
    # [ (name1, ('U', 10)), (name2, ('U', 20)) ]

    for value in dlist.unique_values: # ('U', 10)
        a = Time.now()
        names = dlist.unique_names_for(value)

        datas = np.array([get_fits_data(d, fast=fast) for d in names ])

        if normalize:
            datas = np.array([ d/np.mean(d) for d in datas ])

        if method is 'average':
            combined = np.average(datas, axis=0)
        else:
            combined = np.median(datas, axis=0)
        del datas # saving memory

        print(f'{keys} {value} -> {len(names)} elements.')
        print(f'Done in {Time.now().unix - a.unix :.1f}s')


def correct(keys, pattern, prod, master, operation=None, fast=True):
    '''
    Take a pattern of file names. Correct data against a master.
    Use To subtract or divide.
    '''

    print('Getting headers of all files in pattern')
    dlist = dfits(pattern).fitsort(keys)
    # [ (name1, ('U', 10)), (name2, ('U', 20)) ]

    for value in dlist.unique_values: # ('U', 10)
        a = Time.now()
        names = dlist.unique_names_for(value)

        for name in names:
            if operation != 'flat':
                get_fits_data(name, fast=fast) - master
            else:
                get_fits_data(name, fast=fast) / master

        print(f'{keys} {value} -> {len(names)} elements.')
        print(f'Done in {Time.now().unix - a.unix :.1f}s')


def subtract(keys, pattern, prod, master=None, fast=True):
    '''
    Subtract two images
    '''
    if master is None:
        master = 0
    correct(keys, pattern, prod, master, fast=fast)


def divide(keys, pattern, prod, master=None, fast=True):
    '''
    Divide two images
    '''
    if master is None:
        master = 1
    correct(keys, pattern, prod, master, operation='flat', fast=fast )


def ccdproc_mbias(pattern, output_file=None, header=None, method='median'):
    '''
    CCDproc-based master bias routine.
    Calculates the master bias of a list of fits data.
    Default combining method is median.
    No output file is provided by default.
    '''

    joined_fits = [get_fits_data(p) for p in pattern]
    ccd_data_list = [ccdp.CCDData(d, unit="adu") for d in joined_fits ]

    combined_data = ccdp.combine(ccd_data_list, method=method, unit=u.adu,
                                 dtype=np.uint16, mem_limit=512e6, #)
                                 sigma_clip=True)  #, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 #sigma_clip_func=np.ma.median) #, signma_clip_dev_func=mad_std)
                                 #combined_bias.meta['combined'] = True

    if output_file:
        combined_bias.write(output_file, header=header, overwrite=True)

    return combined_data


def update_keyword(header, key, *tup, comment=None):
    '''
    By Anna Marini
    Updates or create a keyword/value header pair of a given fits file list.
    '''
    value = tup[0].upper()
    time = Time.now()
    hist = time.isot[:-4]+" "

    if key not in header or not header[key]:
        hist += "Created "+key+". "
    else:
        hist += "Updated "+key+". Old value: "+header[key]+". "

    if comment is not None:
        hist += comment

    header[key] = tup
    header.add_history(hist)

    return header


def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # argv[0] is the filename.
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
