#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
from astropy import log
from astropy.io import fits, ascii
from astropy.stats import sigma_clip
from astropy.table import Table
from astropy.time import Time

from astropy.wcs import WCS
from astroquery.mast import Catalogs
from photutils import SkyCircularAperture, SkyCircularAnnulus, aperture_photometry
from astropy.coordinates import SkyCoord

import astropy.units as u
import fitsio
import numpy as np

import cv2

import sorters as s # apparently, no cross imports
from naming import output_file, hist
from fill_header import init_observatory

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

    return 0 if not finfo_list else 1 # finfo=0 # finfo_list[0]=1


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


def master_bias(pattern, keys=[]):
    generic(pattern, keys=keys, method="median", product="MBIAS")

def master_dark(pattern, keys=[], mbias=None):
    generic(pattern, keys=keys, method="median", product="MDARK",
            mbias=mbias)

def master_flat(pattern, keys=[], mbias=None, mdark=None):
    generic(pattern, keys=keys, method="median", product="MFLAT",
            mbias=mbias, mdark=mdark, normalize=True)

def correct_image(pattern, keys=[], mbias=None, mdark=None, mflat=None,
                  method='slice', new_header=False):
    generic(pattern, keys=keys, method=method, product="CLEAN",
            mbias=mbias, mdark=mdark, mflat=mflat, new_header=new_header)

def generic(pattern, keys=[], normalize=False, method=None,
            mbias=None, mdark=None, mflat=None, product=None,
            new_header=False):

    log.info(f'fitsort {len(pattern)} filenames per {keys}')

    df = s.dfits(pattern)
    sortlist = df.fitsort(keys)
    heads = df.heads

    if new_header: o = init_observatory(new_header)

    for value in sortlist.unique_values :
        filenames = sortlist.unique_names_for(value)
        log.info(f'getting {len(filenames)} filenames for {value}')

        # Combine (and save) data per data.
        if method == "slice" or method == "individual":
            for i,filename in enumerate(sorted(filenames)):
                datas = get_fits_data(filename)
                output = combine(filename, normalize=normalize,
                                 mbias=mbias, mdark=mdark, mflat=mflat)

                header = o.newhead(heads[i]) if new_header else heads[i]

                closing(keys, value, product, output, counter=i,
                        header=header)

        # Combine and save acting on a data cube
        else:
            datas = np.array([ get_fits_data(f) for f in filenames ])
            output = combine(datas, normalize=normalize, method=method,
                             mbias=mbias, mdark=mdark, mflat=mflat)

            header = o.newhead(header=heads[0]) if new_header else heads[0]

            closing(keys, value, product, output, header=header)


def closing(keys, value, product, output, counter=False, header=False):

    if header:
        #header = heads[0].copy() # TODO choose head per head
        header.add_history(hist())

    text = dict(zip(keys, value)) if keys else None
    outfile =  output_file(product=product, text=text, counter=counter)
    write_fits(output, outfile, header=header, fast=False)


def combine(datas, normalize=False, method=None, precision='float32',
            mbias=None, mdark=None, mflat=None, mask=False):

    a = Time.now()

    # Datas from pattern
    if isinstance(datas, str): datas = [datas]
    if isinstance(datas, list) and isinstance(datas[0], str):
        datas = np.array([get_fits_data(d) for d in datas ])

    # Master datas from filename
    if isinstance(mbias, str): mbias = get_fits_data(mbias)
    if isinstance(mdark, str): mdark = get_fits_data(mdark)
    if isinstance(mflat, str): mflat = get_fits_data(mflat)

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

    if method == 'average':
        combined = np.average(datas, axis=0).astype(precision)
    elif method == 'median':
        combined = np.median(datas, axis=0).astype(precision)
    else: # cube or 1-slice cube.
        combined = np.squeeze(datas)

    log.info(f'{method}: {datas.shape}{datas.dtype} -> {combined.shape}{combined.dtype}')
    log.info(f'Done in {Time.now().unix - a.unix :.1f}s')
    del datas # Saving memory

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


def detect_sources(image):
    ''' By Anna Marini
    Extract the light sources from the image
    '''
    threshold = detect_threshold(image, nsigma=2.)
    sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    mean, median, std = sigma_clipped_stats(image, sigma=3)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(image - median)
   # Pixel coordinates of the sources
    x = np.array(sources['xcentroid'])
    y = np.array(sources['ycentroid'])
    return x,y


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
    '''
    Use opencv to find centroids of highly defocused images template matching.
    '''

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


def load_catalog(filename=False, header=False, wcs=False, ra_key=False, dec_key=False):
    '''
    From Anna Marini: get positions from catalog.
    '''

    if filename and not header: header = get_fits_header(filename)
    if header and not wcs:
        wcs = WCS(header)
        if ra_key and dec_key:
            ra = header[ra_key]
            dec = header[dec_key]

    ra = wcs.wcs.crval[0]
    dec = wcs.wcs.crval[1]

    # Diagonal
    diag_bound = wcs.pixel_to_world_values([ [0,0], wcs.pixel_shape ])
    radius = np.mean(diag_bound[1] - diag_bound[0]) / 2

    catalog = Catalogs.query_region(f'{ra} {dec}',
                                    frame='fk5',
                                    unit="deg",
                                    radius=radius,
                                    catalog = 'TIC')

    return catalog


def set_apertures(catalog, limit=16, r=10, r_in=18, r_out=22):
    '''From Anna Marini: get a catalog and
    set apertures and annulus for photometry.
    '''
    radec = catalog['ra','dec','Vmag']
    mask = radec['Vmag']<limit
    radec = radec[mask]

    positions = SkyCoord(radec['ra'], radec['dec'],
                         frame='fk5',
                         unit=(u.deg,u.deg))

    aperture = SkyCircularAperture(positions,
                                   r=r*u.arcsec)
    annulus = SkyCircularAnnulus(positions,
                                 r_in=r_in*u.arcsec,
                                 r_out=r_out*u.arcsec)
    apers = [aperture, annulus]

    return apers


def do_photometry(data, apers, wcs, obstime=False):

    phot_table = aperture_photometry(data, apers, wcs=wcs)

    pixar = apers[0].to_pixel(wcs)
    pixan = apers[1].to_pixel(wcs)

    bkg_mean = phot_table['aperture_sum_1'] / pixan.area
    bkg_sum = bkg_mean * pixar.area
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    
    phot_table['residual_aperture_sum'] = final_sum
    phot_table['mjd-obs'] = obstime
    
    return(phot_table)


def apphot(pattern, display=False):

    pattern = sorted(pattern)

    header0 = get_fits_header(pattern[0])
    wcs0 = WCS(header0)
    
    catalog = load_catalog(wcs=wcs0)
    apers = set_apertures(catalog)

    tables = Table()

    if display:
        import pyds9
        d=pyds9.DS9("ds9")
        
    for filename in pattern:
        header = get_fits_header(filename)
        data = get_fits_data(filename)
        wcs = WCS(header)
    
        #catalog = load_catalog(wcs=wcs)
        #apers = set_apertures(catalog)
            
        phot_table = do_photometry(data, apers, wcs, obstime=header["MJD-OBS"])

        positions = SkyCoord(catalog['ra'], catalog['dec'],
                             frame='fk5',
                             unit=(u.deg,u.deg))                

        if display:
            d.set(f"file {filename}")
            
            # for i,pos in enumerate(positions):
            #     p = pos.to_pixel(wcs)
            #     circ = f'circle({p[0]}, {p[1]}, {10})'
            #     d.set("regions", circ)
            #     d.set("region", f"text {p[0]} {p[1]} "+"{"+str(i)+"}")
                
            for i,a in enumerate(apers[0].to_pixel(wcs)) :
                circ = f'circle({a.positions[0]}, {a.positions[1]}, {a.r})'
                d.set("regions", circ)
                d.set("region", f"text {a.positions[0]}, {a.positions[1]} "+"{"+str(i)+"}")
            
            for a in apers[1].to_pixel(wcs):
                circ = f'circle({a.positions[0]}, {a.positions[1]}, {a.r_in})'
                d.set("regions", circ)
                circ = f'circle({a.positions[0]}, {a.positions[1]}, {a.r_out})'
                d.set("regions", circ)
        
        tables.add_column(phot_table["residual_aperture_sum"], rename_duplicate=True)
            
        log.info(f"Done {filename}")

    return tables


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
