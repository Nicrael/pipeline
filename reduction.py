#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy import log
from astropy.io import fits, ascii
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
import fitsio
import numpy as np

import cv2

import sorters as s # apparently, no cross imports
from naming import output_file, hist

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

    log.debug(f"Getting header from {filename}")
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

    log.debug(f"Getting data from {filename}")
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


def write_fits(data, output_file, header=None, fast=False):
    '''
    Write a fits file.
    It adds a checksum keyword.
    '''

    if fast:
        hdu = fitsio.FITS(output_file,'rw')
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


def master_bias(pattern, keys):
    generic(pattern, keys=keys, method="median", product="MBIAS")

def master_flat(pattern, keys, mbias):
    generic(pattern, keys=keys, method="median", product="MFLAT",
            mbias=mbias, normalize=True)

def correct_image(pattern, keys, mbias, mflat, method='slice'):
    generic(pattern, keys, method=method, product="CLEAN",
            mbias=mbias, mflat=mflat)

def generic(pattern, keys=[], normalize=False, method=None,
            mbias=None, mdark=None, mflat=None, product=None):

    log.info(f'fitsort {len(pattern)} files per {keys}')

    df = s.dfits(pattern)
    sortlist = df.fitsort(keys)
    heads = df.heads

    for value in sortlist.unique_values :
        files = sortlist.unique_names_for(value)
        log.info(f'getting {len(files)} files for {value}')

        # Combine (and save) data per data.
        if method is "slice" or method is "individual":
            for i,f in enumerate(files):
                datas = get_fits_data(f)
                output = combine(f, normalize=normalize,
                                 mbias=mbias, mdark=mdark, mflat=mflat)

                closing(heads, keys, value, product, output, counter=i)

        # Combine and save acting on a data cube
        else:
            datas = np.array([ get_fits_data(f) for f in files ])
            output = combine(datas, normalize=normalize, method=method,
                             mbias=mbias, mdark=mdark, mflat=mflat)

            closing(heads, keys, value, product, output)


def closing(heads, keys, value, product, output, counter=False):

    header = heads[0].copy() # TODO choose head per head
    header.add_history(hist())
    text = dict(zip(keys, value)) if keys else None
    outfile =  output_file(product=product, text=text, counter=counter)
    write_fits(output, outfile, header=header, fast=False)


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

        # if normalize:
        #     for i, d in enumerate(datas):
        #         datas[i] = (d/np.mean(d)).astype(precision)

        if method is 'average':
            combined = np.average(datas, axis=0).astype(precision)
        elif method is 'median':
            combined = np.median(datas, axis=0).astype(precision)
        else: # cube or 1-slice cube.
            combined = np.squeeze(datas)

        log.info(f'{method}: {datas.shape}{datas.dtype} -> {combined.shape}{combined.dtype}')

        del datas # Saving memory

    else: # len(datas) == 0
        log.warn('Input datas are empty:  Result is input.')
        combined = np.array(datas, dtype=np.uint16)

    log.info(f'Done in {Time.now().unix - a.unix :.1f}s')
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



def detect_sources(pattern):
    ''' By Anna Marini
    Extract the light sources from the image
    '''
    rot = rotate_img(pattern)
    threshold = detect_threshold(rot, nsigma=2.)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    mean, median, std = sigma_clipped_stats(rot, sigma=3)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(rot - median)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
   # Pixel coordinates of the sources
    d = dict();
    d['x1'] = np.array(sources['xcentroid'])
    d['y1'] = np.array(sources['ycentroid'])
    return(d)


def rescale(array):
    '''
    Takes an array.
    Rescale min to 0,  max to 255, then
    change dtype, as opencv loves uint8 data type.
    Returns the rescaled uint8 array.
    '''
    array -= np.min(array)
    array = array/(np.max(array)/255.0)
    return array.astype(np.uint8)


def detect_donuts(filename, template):

    img = rescale(get_fits_data(filename))
    tpl = rescale(get_fits_data(template))

    res = cv2.matchTemplate(img, tpl, cv2.TM_CCOEFF_NORMED)
    threshold = 0.6

    loc = np.where( res >= threshold)
    x,y = loc
    p = np.repeat("point ", y.size)
    t = [p, (y+tpl.shape[0]/2),(x+tpl.shape[1]/2)]
    table = Table(t, names=['# ','## ','###']) # bleah
    ascii.write(table, "donuts.reg", overwrite=True)

    return res


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
