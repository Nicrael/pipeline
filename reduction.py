#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy.io import fits, ascii
from astropy.time import Time
import astropy.units as u
import numpy as np
from astropy.stats import sigma_clip

from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std
from astropy.table import Table

#                  1               2            3                   4                 5           6
# biases ------> MBIAS
#        darks - MBIAS ->   darks_debiased
#                           darks_debiased -> MDARK
#        flats - MBIAS ->   flats_debiased
#                           flats_debiased  - MDARK -->  flats_debiased_dedarked
#                                                        flats_debiased_dedarked -> MFLAT
#      objects - MBIAS -> objects_debiased
#                         objects_debiased  - MDARK -> objects_debiased_dedarked
#                                                      objects_debiased_dedarked  / MFLAT -> objects_reduc
#

def choose_hdu(filename):
    '''
    Detect wether the fits file is compressed with
    fpack, and choose the right HDU.
    '''
    finfo = fits.info(filename, output=False) # Returns a list of tuples.
    finfo_list = [item for item in finfo if 'COMPRESSED_IMAGE' in item]
    if not finfo_list:
        return 0 # finfo[0][0] # 0 if not compressed
    else:
        return 1 # finfo_list[0][0] # 1 if compressed


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


def frame_dict(filename, with_data=False):
    '''
    Create a dictionary related to an observation frame
    '''
    fd = {
        'name' : filename,
        'head' : get_fits_header(filename),
        'data' : None,
        }
    if with_data:
        fd['data'] = get_fits_data(filename)

    return fd


class AttrDict(dict):
    '''
    Create an objects where properties are dict keys.
    '''
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def frame(filename):
    '''
    Create a frame object related to an observation frame
    '''
    fr = AttrDict(frame_dict(filename))

    return fr


def frame_list(pattern):
    '''
    Create a list of frames from filename pattern
    '''
    list1 = []
    for filename in pattern:
        list1.append(frame(filename))

    return list1


def join_fits_header(pattern):
    '''
    Join the header of list of fits files in a list.
    '''
    heads = np.array([ get_fits_header(f) for f in pattern ])
    return heads


def join_fits_data(pattern):
    '''
    Join the data of a list of fits files in a tuple.
    Tuple format is useful for stacking in a data cube.
    '''
    datas = np.array([ get_fits_data(f) for f in pattern ])
    return datas


def stack_fits_data(datas):
    '''
    Stack a list of fits datas in a data cube.
    It is useful to perform pixel-per-pixel operations,
    such as an average.
    '''
    datacube = np.dstack(datas)
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
    datatype=datacube.dtype
    average = np.average(datacube, axis=2)
    return average.astype(datatype)


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


def oarpaf_combine(pattern, method='median', output_file=None, header=None):
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

    if output_file:
        write_fits(combined_data, output_file, header=header)

    return combined_data


def oarpaf_mask(data, sigma=3, output_file=None, header=None):

    mask = sigma_clip(data, masked=True).mask.astype(int)

    if output_file:
        write_fits(mask, output_file, header=header)

    return mask


def oarpaf_mask_reg(data, sigma=3, output_file=None):

    y,x = np.where(data == True)
    p = np.repeat("point ", y.size)
    t = [p, x+1,y+1]

    table = Table(t, names=['# ','## ','###']) # bleah

    if output_file:
        ascii.write(table, output_file, overwrite=True)

    return table


def ccdproc_mbias(pattern, output_file=None, header=None, method='median'):
    '''
    CCDproc-based master bias routine.
    Calculates the master bias of a list of fits data.
    Default combining method is median.
    No output file is provided by default.
    '''

    joined_fits = join_fits_data(pattern)
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

    if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()

#### Cemetery of old functions #####

# def new_header():
#     return fits.PrimaryHDU().header

# def to_list(arg):
#     if type(arg) is not list: arg = [ arg ]
#     return arg

# def is_number(s):
#     '''
#     Check if a string contains a (float) number.
#     Useful to test decimal or sexagesimal coordinates.
#     '''
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False


# def to_number(s):
#     try:
#         return int(s)
#     except ValueError:
#         return float(s)


# def get_fits_data_or_header(filename,get):
#     '''
#     Return the header or the data of a fits file.
#     '''
#     which_hdu = choose_hdu(filename)
#     with fits.open(filename) as hdul:
#         if get is 'header':
#             return hdul[which_hdu].header;
#         elif get is 'data':
#             return hdul[which_hdu].data;
#         else:
#             return


# def get_fits_data2(filename):
#     '''
#     Return the data of the fits file.
#     Alternative method based on fitsio.
#     '''
#     which_hdu = choose_hdu(filename)
#     with fitsio.FITS(filename) as f:
#         data = f[which_hdu].read()
#         return data


# From nested dict (json) to object
# class obj(object):
#     def __init__(self, d):
#         for a, b in d.items():
#             if isinstance(b, (list, tuple)):
#                 setattr(self, a, [obj(x) if isinstance(x, dict) else x for x in b])
#             else:
#                 setattr(self, a, obj(b) if isinstance(b, dict) else b)
