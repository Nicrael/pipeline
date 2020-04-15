#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.rotations import rotation_matrix
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np


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
    Join the header of a given fits file list in a tuple.
    '''
    tuple_of_fits_header = ()
    for filename in pattern:
        fits_header = get_fits_header(filename)
        tuple_of_fits_header += (fits_header,)
    return tuple_of_fits_header


def join_fits_data(pattern):
    '''
    Join the data of a given fits file list in a tuple.
    It is a useful format for stacking in a data cube.
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


def oarpaf_mbias(pattern, output_file=None, method='median'):
    '''
    Custom master bias routine.
    Calculates the master bias of a list of of fits files.
    Default combining method is median.
    No output file is provided by default.
    '''
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

    

# def ccdproc_mbias(pattern, output_file=None, method='median'):
#     '''
#     CCDproc-based master bias routine.
#     Calculates the master bias of a list of of fits files.
#     Default combining method is median.
#     No output file is provided by default.
#     '''
#     combined_data = ccdp.combine(pattern, method=method, unit=u.adu,
#                                  dtype=np.uint16, mem_limit=150e6)
#     # sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
#     # sigma_clip_func=np.ma.median, signma_clip_dev_func=mad_std)
#     #combined_bias.meta['combined'] = True
#     if output_file is not None:
#         combined_bias.write(output_file, overwrite=True)
#     return combined_data


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
    header['history'] = hist

    return header


def init_wcs(header):
    '''
    By Anna Marini
    Valid for SPM 84cm FITS headers.
    Updates header with basic WCS keywords to help astrometry.net
    '''
    angle = 90 * u.deg
    plate = 0.25 * header['CCDXBIN']*u.arcsec

    rot_matrix = plate.to(u.deg) * rotation_matrix(angle)[:-1,:-1]
    coor = SkyCoord(header['RA'], header['DEC'],
                    unit=(u.hourangle, u.deg),
                    obstime=Time(header['JD'], format='jd'))
    w = WCS(header)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [header['NAXIS1']/2, header['NAXIS2']/2]
    w.wcs.crval = [coor.ra.degree, coor.dec.degree]
    w.wcs.cd = rot_matrix.to(u.deg).value

    header.extend(w.to_header(), update=True)    

    return header

    
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
