#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.time import Time
import astropy.units as u
import numpy as np

from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std


#                  1               2            3                   4                 5           6
# biases ------> MBIAS
#        darks - MBIAS ->   darks_debiased
#                           darks_debiased -> MDARK
#        flats - MBIAS ->   flats_debiased
#                           flats_debiased  - MDARK -->  flats_debiased_dedarked -> MFLAT
#      objects - MBIAS -> objects_debiased
#                         objects_debiased  - MDARK -> objects_debiased_dedarked  / MFLAT -> objects_reduc
#

# biases, darks, flat
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 3 MDARK               = combine(darks_debiased)
# 4 *_debiased_dedarked = subtract_dark(*_debiased, MDARK)
# 5 MFLAT               = combine(flats_debiased_dedarked)
# 6 objects_reduc       = flat_correct(objects_debiased_dedarked, MFLAT)
#

# biases,      , flat
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 5 MFLAT               = combine(flats_debiased)
# 6 objects_reduc       = flat_correct(objects_debiased, MFLAT)
#

# biases, darks
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 3 MDARK               = combine(darks_debiased)
# 4 *_debiased_dedarked = subtract_dark(*_debiased, MDARK)
#

#         darks, flat
#
# 3 MDARK               = combine(darks)
# 4 *_dedarked          = subtract_dark(*, MDARK)
# 5 MFLAT               = combine(flats_dedarked)
# 6 objects_reduc       = flat_correct(objects_dedarked, MFLAT)
#

# biases
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
#

#         darks
#
# 3 MDARK               = combine(darks)
# 4 *_dedarked          = subtract_dark(*, MDARK)
#

#                 flat
#
# 5 MFLAT               = combine(flats)
# 6 objects_reduc       = flat_correct(objects, MFLAT)
#


def choose_hdu(filename):
    '''
    Detect wether the fits file is compressed with
    fpack, and choose the right HDU.
    '''
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return finfo[0][0] # 0 if not compressed
    else:
        return finfo_list[0][0] # 1 if compressed


def get_fits_data_or_header(filename,get):
    '''
    Return the header or the data of a fits file.
    '''
    which_hdu = choose_hdu(filename)
    with fits.open(filename) as hdul:
        if get is 'header':
            return hdul[which_hdu].header;
        elif get is 'data':
            return hdul[which_hdu].data;
        else:
            return


def get_fits_header(filename):
    '''
    Return the header of the fits file.
    '''
    #return get_fits_data_or_header(filename,'header')
    which_hdu = choose_hdu(filename)
    return fits.getheader(filename, which_hdu)


def get_fits_data(filename):
    '''
    Return the data of the fits file.
    '''
    #return get_fits_data_or_header(filename,'data')
    which_hdu = choose_hdu(filename)
    return fits.getdata(filename, which_hdu)


# def get_fits_data2(filename):
#     '''
#     Return the data of the fits file.
#     Alternative method based on fitsio.
#     '''
#     which_hdu = choose_hdu(filename)
#     with fitsio.FITS(filename) as f:
#         data = f[which_hdu].read()
#         return data


def join_fits_header(pattern):
    '''
    Join the header of list of fits files in a tuple.
    '''
    tuple_of_fits_header = ()
    for filename in pattern:
        fits_header = get_fits_header(filename)
        tuple_of_fits_header += (fits_header,)
    return tuple_of_fits_header


def join_fits_data(pattern):
    '''
    Join the data of a list of fits files in a tuple.
    Tuple format is useful for stacking in a data cube.
    '''
    tuple_of_fits_data = ()
    for filename in pattern:
        fits_data = get_fits_data(filename)
        tuple_of_fits_data += (fits_data,)
    return tuple_of_fits_data


def stack_fits_data(tuple_of_fits_data):
    '''
    Stack the fits data in a data cube.
    It is useful to perform pixel-per-pixel operations,
    such as an average.
    '''
    datacube = np.dstack(tuple_of_fits_data)
    return datacube


def median_datacube(datacube):
    '''
    Make a median of a data cube.
    '''
    datatype=datacube.dtype
    median = np.median(datacube, axis=2, overwrite_input=True)
    return median.astype(datatype)


def average_datacube(datacube):
    '''
    Make an average of a data cube.
    '''
    average = np.average(datacube, axis=2)
    return average


def write_fits(data, filename, header=None):
    '''
    Writes a fits file.
    It adds a checksum keyword.
    '''
    if header is None:
        hdu = fits.PrimaryHDU(data)
    else:
        hdu = fits.PrimaryHDU(data, header=header)
    hdu.writeto(filename, overwrite=True, checksum=True)
    return hdu


def oarpaf_mbias(pattern, method='median', output_file=None, header=None):
    '''
    Custom master bias routine.
    Calculates the master bias of a list of of fits files.
    Default combining method is median.
    No output file is provided by default.
    '''
    joined_fits = join_fits_data(pattern)
    datacube = stack_fits_data(joined_fits)
    del joined_fits # saving memory
    if method is 'average':
        combined_data = average_datacube(datacube)
    else:
        combined_data = median_datacube(datacube)
    del datacube # saving memory

    header = get_fits_header(pattern[0]) if header else None
    if output_file is not None:
        write_fits(combined_data, output_file, header=header)
    return combined_data


def ccdproc_mbias(pattern, output_file=None, method='median'):
    ccd_data_list = [ccdp.CCDData(p, unit="adu")  for p in pattern ]
    '''
    CCDproc-based master bias routine.
    Calculates the master bias of a tuple of fits data.
    Default combining method is median.
    No output file is provided by default.
    '''

    combined_data = ccdp.combine(ccd_data_list, method=method, unit=u.adu,
                                 dtype=np.uint16, mem_limit=1024e6)
    # sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
    # sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std)
    #combined_bias.meta['combined'] = True
    if output_file is not None:
        combined_bias.write(output_file, overwrite=True)
    return combined_data


def is_number(s):
    '''
    Check if a string contains a (float) number.
    Useful to test decimal or sexagesimal coordinates.
    '''
    try:
        float(s)
        return True
    except ValueError:
        return False


def to_number(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def xyval(arr):
    '''
    Transform a matrix of values in a "x,y,value table".
    Slow.
    '''
    y,x = np.indices(arr.shape)
    return x.ravel()+1, y.ravel()+1, arr.ravel()


def xyfilter(arr,value=True):
    '''
    Filter where the value of an "x,y,value table" is equal to a specific value.
    Default value is True. Useful to calculate a bad pixel map for ds9.
    '''
    return arr[0][arr[2]],arr[1][arr[2]]


# def update_keyword(key,valore,pattern):
#     '''
#     By Anna Marini
#     Updates or create a keyword/value header pair of a given fits file list.
#     '''
#     for filename in pattern:
#         which_hdu = choose_hdu(pattern)
#         with fits.open(filename,'update') as hdul:
#             header = hdul[which_hdu].header;
#             header[key] = valore
#             return header

# From nested dict (json) to object
# class obj(object):
#     def __init__(self, d):
#         for a, b in d.items():
#             if isinstance(b, (list, tuple)):
#                 setattr(self, a, [obj(x) if isinstance(x, dict) else x for x in b])
#             else:
#                 setattr(self, a, obj(b) if isinstance(b, dict) else b)



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


def new_header():
    return fits.PrimaryHDU().header

def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".

    # import datetime

    # In [357]: clipped=stats.sigma_clip(ccdbias,sigma=3)
    # In [358]: bbb=mbias.xyval(clipped.mask)
    # In [359]: ddd=mbias.xyfilter(bbb)
    # In [360]: ascii.write(ddd,'mask.reg',overwrite=True)

    # dfits ~/desktop/oarpaf/test-sbig-stx/*.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'

    #a = datetime.datetime.now()
    oarpaf_mbias(pattern)
    #b = datetime.datetime.now()
    #c = b - a ; print(c.total_seconds())

    # a = datetime.datetime.now()
    # ccdproc_mbias(pattern)
    # b = datetime.datetime.now()
    # c = b - a ; print(c.total_seconds())


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()


### SBAGLIATO ###

#                  1            2                  3          4,5         6        7
# biases ------> MBIAS
#        darks - MBIAS -> darks_clean
#                         darks_clean --------> MDARK
#                                       flats - MDARK -> flats_clean --> MFLAT
#                                     objects - MDARK -> objects_clean / MFLAT -> objects_reduc

# biases, darks, flat
#
# 1 MBIAS         = combine(biases)
# 2 darks_clean   = subtract_bias(darks, MBIAS)
# 3 MDARK         = combine(darks_clean)
# 4 flats_clean   = subtract_dark(flats, MDARK)
# 5 objects_clean = subtract_dark(objects, MDARK)
# 6 MFLAT         = combine(flats_clean)
# 7 objects_reduc = flat_correct(objects_clean, MFLAT)

# biases, darks
#
# 1 MBIAS         = combine(biases)
# 2 darks_clean   = subtract_bias(darks, MBIAS)
# 3 MDARK         = combine(darks_clean)
# 5 objects_clean = subtract_dark(objects, MDARK)
#

# biases,        flat
#
# 1 MBIAS         = combine(biases)
# 4 flats_clean   = subtract_dark(flats, MBIAS)
# 5 objects_clean = subtract_dark(objects, MBIAS)
# 6 MFLAT         = combine(flats_clean)
# 7 objects_reduc = flat_correct(objects_clean, MFLAT)
#

#         darks, flat
#
# 3 MDARK         = combine(darks)
# 4 flats_clean   = subtract_dark(flats, MDARK)
# 5 objects_clean = subtract_dark(objects, MDARK)
# 6 MFLAT         = combine(flats_clean)
# 7 objects_reduc = flat_correct(objects_clean, MFLAT)
#

# biases
#
# 1 MBIAS         = combine(biases)
# 5 objects_clean = subtract_dark(objects, MBIAS)
#

#         darks, flat
#
# 3 MDARK         = combine(darks)
# 5 objects_clean = subtract_dark(objects, MDARK)
#

#                flat
#
# 6 MFLAT         = combine(flats)
# 7 objects_reduc = flat_correct(objects_clean, MFLAT)
#
