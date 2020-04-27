#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules

from astropy.io import fits, ascii
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.stats import sigma_clip
from astropy.table import Table, unique
from astropy.time import Time
import astropy.units as u
import ccdproc as ccdp
import fitsio
import numpy as np

##########################################################################

#heads = [ get_fits_header(f, fast=True) for f in pattern ]
#table = [ dict([ [h["name"],h["value"]] for h in H.records()]) for H in heads ]

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

    def unique_for(self, keys):
        # if isinstance(keys, str):
        #     keys = [keys]
        self.data = self.table.group_by(keys)
        self.unique = self.data.groups.keys.as_array().tolist()
        return self.unique

    def names_for(self, keys):
        if isinstance(keys, str):
            keys = [keys]
        self.names = [ np.array(g[keys]).tolist() for g in self.data.groups]
        self.data = self.table[keys]
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


def get_fits_data(filename, fast=True):
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


def load(pattern):
    '''
    Return the datas of a file pattern.
    '''
    if isinstance(pattern, str): pattern = [pattern]
    datas = np.array([ get_fits_data(f) for f in pattern ])

    # Collapsing array of cubes (3,30,100,100) -> (90,100,100)
    if len(datas.shape) > 3 :
        datas = datas.reshape(-1, *datas.shape[-2:])

    return datas.squeeze()


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


def write_fits(data, output_file, header=None):
    '''
    Write a fits file.
    It adds a checksum keyword.
    '''

    if header is None:
        hdu = fits.PrimaryHDU(data)
    else:
        hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(output_file, overwrite=True, checksum=True)
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



def generic(pattern, keys=[], normalize=False, method=None,
            mbias=None, mdark=None, mflat=None, output_file=None):

    print(f'fitsort {len(pattern)} files per {keys}')

    sortlist = minidb(pattern).group_by(keys)

    for value in sortlist.unique :
        files = sortlist.names_for(value)
        print(f'getting {len(files)} files for {value}')

        # Combine (and save) data per data.
        if method is "slice" or method is "individual":
            for i,f in enumerate(files):
                datas = get_fits_data(f)
                output = combine(f, normalize=normalize,
                                 mbias=mbias, mdark=mdark, mflat=mflat)

                outfile =  f'generic{value}{i}.fits' if not output_file else output_file
                print(outfile)
                write_fits(output, outfile)

        # Combine and save acting on a data cube
        else:
            datas = np.array([get_fits_data(f) for f in files ])
            output = combine(datas, normalize=normalize, method=method,
                             mbias=mbias, mdark=mdark, mflat=mflat)

            outfile =  f'generic{value}.fits' if not output_file else output_file
            print(outfile)
            write_fits(output, outfile)


def combine(datas, normalize=False, method=None,
            mbias=None, mdark=None, mflat=None, mask=False):

    a = Time.now()

    if len(datas): # Check if datas is not empty

        # Datas from pattern
        if isinstance(datas, str): datas = [datas]
        if isinstance(datas, list) and isinstance(datas[0], str):
            datas = np.array([get_fits_data(d) for d in datas ])

        # Master datas from filename
        if isinstance(mbias, str): mbias = get_fits_data(mbias)
        if isinstance(mdark, str): mdark = get_fits_data(mdark)
        if isinstance(mflat, str): mflat = get_fits_data(mflat)

        precision='float32'

        # Cannot cast type
        if mbias is not None and len(mbias):
            datas = (datas - mbias).astype(precision)
        if mdark is not None and len(mdark):
            datas = (datas - mdark).astype(precision)
        if mflat is not None and len(mflat):
            datas = (datas / mflat).astype(precision)

        del mbias, mdark, mflat

        # Did not find a faster method to save memory.
        if normalize:
            bottle = np.zeros(shape=datas.shape).astype(precision)
            for i, d in enumerate(datas):
                bottle[i] = d/np.mean(d).astype(precision)
            datas = bottle
            del bottle

        if method is 'average':
            combined = np.average(datas, axis=0).astype(precision)
        elif method is 'median':
            combined = np.median(datas, axis=0).astype(precision)
        else: # cube or 1-slice cube.
            combined = np.squeeze(datas)

        print(f'{method}: {datas.shape}{datas.dtype} -> {combined.shape}{combined.dtype}')

        del datas # Saving memory

    else: # len(datas) == 0
        print('Input datas are empty:  Result is input.')
        combined = np.array(datas, dtype=np.uint16)

    print(f'Done in {Time.now().unix - a.unix :.1f}s')
    return combined


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
