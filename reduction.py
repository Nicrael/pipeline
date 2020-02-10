#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.io import fits
import fitsio
import datetime

from astropy.nddata import CCDData
from astropy.stats import mad_std
import astropy.units as u

from pathlib import Path
import ccdproc as ccdp
import matplotlib.pyplot as plt

'''
Detect is the fits file is compressed with fpack, and choose the right HDU.
'''
def choose_hdu(filename):
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0] # 1 if compressed

'''
Return the header of the fits file.
'''
def get_fits_header(filename):
    which_hdu = choose_hdu(filename)
    with fits.open(filename) as hdul:
        header = hdul[which_hdu].header;
        return header

'''
Return the data of the fits file.
'''
def get_fits_data(filename):
    which_hdu = choose_hdu(filename)
    with fits.open(filename) as hdul:
        data = hdul[which_hdu].data;
        return data
    # with fitsio.FITS(filename) as f:
    #     data = f[which_hdu].read()
    #     return data
    
'''
Join the data of a given fits file list in a tuple.
It is a useful format for stacking in a data cube
'''
def join_fits_data(pattern):
    tuple_of_fits_data = ()
    for filename in pattern:
        fits_data = get_fits_data(filename)
        tuple_of_fits_data += (fits_data,)
    return tuple_of_fits_data

'''
Stack the fits data in a data cube.
It is useful to perform pixel-per-pixel operations such as an average.
'''
def stack_fits_data(tuple_of_fits_data):
    datacube = np.dstack(tuple_of_fits_data)
    return datacube

'''
Make a median of a data cube.
'''
def median_datacube(datacube):
    datatype=datacube.dtype
    median = np.median(datacube, axis=2, overwrite_input=True)
    return median.astype(datatype)

'''
Make an average of a data cube.
'''
def average_datacube(datacube):
    average = np.average(datacube, axis=2)
    return average

'''
Writes a fits file.
It adds a checksum keyword.
'''
def write_fits(data, filename, header=None):
    if header is None:
        hdu = fits.PrimaryHDU(data)
    else:
        hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(filename, overwrite=True, checksum=True)
    return hdu

'''
Custom master bias routine.
Calculates the master bias of a list of of fits files.
Default combining method is median.
No output file is provided by default
'''
def oarpaf_mbias(pattern, output_file=None, method='median'):
    joined_fits = join_fits_data(pattern)
    datacube=stack_fits_data(joined_fits)
    del joined_fits # saving memory
    if method is 'average':
        combined_data = average_datacube(datacube)
    else:
        combined_data = median_datacube(datacube)
    del datacube # saving memory
    header=get_fits_header(pattern[0])
    if output_file is not None:
        write_fits(combined_data, output_file, header=header)
    return combined_data
    
'''
CCDproc-based master bias routine.
Calculates the master bias of a list of of fits files.
Default combining method is median.
No output file is provided by default.
'''
def ccdproc_mbias(pattern, output_file=None, method='median'):

    combined_data = ccdp.combine(pattern, method=method, unit=u.adu, dtype=np.uint16, mem_limit=150e6)
    # sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
    # sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std)
    #combined_bias.meta['combined'] = True
    if output_file is not None:
        combined_bias.write(output_file, overwrite=True)
    return combined_data


# In [357]: clipped=stats.sigma_clip(ccdbias,sigma=3)
# In [358]: bbb=mbias.xyval(clipped.mask)
# In [359]: ddd=mbias.xyfilter(bbb)
# In [360]: ascii.write(ddd,'mask.reg',overwrite=True)


'''
Transform a matrix of values in a "x,y,value table". 
'''
def xyval(arr):
    y,x = np.indices(arr.shape)
    return x.ravel()+1, y.ravel()+1, arr.ravel()
    
'''
Filter where the value of an "x,y,value table" is equal to a specific value.
Default value is True. Useful to calculate a bad pixel map for ds9.
'''
def xyfilter(arr,value=True):
    return arr[0][arr[2]],arr[1][arr[2]]


def main():
    pattern = sys.argv[1:]    # File(s). "1:" significa "dal 3 in poi".

    # dfits *.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'

    a = datetime.datetime.now()
    oarpaf_mbias(pattern)
    b = datetime.datetime.now()
    c = b - a ; print(c.total_seconds())

    # a = datetime.datetime.now()
    # ccdproc_mbias(pattern)
    # b = datetime.datetime.now()
    # c = b - a ; print(c.total_seconds())

if __name__ == '__main__':
    import sys

    # if len(sys.argv) < 4 :    # C'è anche lo [0] che è il nome del file :)
    #     print("Usage:  "+sys.argv[0]+" KEYWORD value <list of FITS files>")
    #     sys.exit()

    main()
