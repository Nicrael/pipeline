#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
FITS file related functions.
'''

# System modules
from astropy import log
from astropy.io import fits

try:
    import fitsio
    FAST = True
except ImportError:
    log.warning("fitsio module not found: cannot use fast mode.")
    FAST = False

# Local modules


def choose_hdu(filename, fast=False):
    '''
    Detect whether the fits file is compressed with
    fpack, and choose the right HDU.
    fast: Alternative mode based on fitsio
    '''
    if fast:
        finfo = fitsio.FITS(filename)  # Object
        finfo_list = [f.get_extnum() for f in finfo if f.is_compressed()]
    else:
        finfo = fits.info(filename, output=False)  # List of tuples.
        finfo_list = [f[0] for f in finfo if 'COMPRESSED_IMAGE' in f]

    return 0 if not finfo_list else 1  # finfo=0 # finfo_list[0]=1


def get_fits_header(filename, fast=False):
    '''
    Return the header of the fits file.
    fast uses fitsio.
    '''
    which_hdu = choose_hdu(filename, fast=fast)
    if fast:
        #header = fitsio.read_header(filename, which_hdu)
        with fitsio.FITS(filename) as file_:
            header = file_[which_hdu].read_header()
    else:
        header = fits.getheader(filename, which_hdu)

    log.debug("Getting header from {filename}", filename=filename)
    return header


def get_fits_data(filename, fast=FAST):
    '''
    Return the data of the fits file.
    If fitsio=True, use fitsio.
    '''
    which_hdu = choose_hdu(filename, fast=fast)
    if fast:
        #data = fitsio.read(filename, which_hdu)
        with fitsio.FITS(filename) as file_:
            data = file_[which_hdu].read()
    else:
        data = fits.getdata(filename, which_hdu)

    log.debug("Getting data from {filename}", filename=filename)
    return data


def write_fits(data, output_file, header=None, fast=False):
    '''
    Write a fits file.
    It adds a checksum keyword.
    '''

    if fast:
        hdu = fitsio.FITS(output_file, 'rw')
        if header:
            hdu.write(data=data, header=header)
        else:
            hdu.write(data=data)
        hdu.close()
    else:
        if header:
            hdu = fits.PrimaryHDU(data, header=header)
        else:
            hdu = fits.PrimaryHDU(data)
        hdu.writeto(output_file, overwrite=True, checksum=True)
        
    log.info(f"Writing fits file to {output_file}")
    return hdu
