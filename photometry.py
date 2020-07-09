#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# System modules
import json
from astropy import log
from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.mast import Catalogs
from photutils import SkyCircularAperture, SkyCircularAnnulus, aperture_photometry
from photutils import DAOStarFinder
from photutils import make_source_mask
import astropy.units as u
import cv2
import numpy as np
import matplotlib.pyplot as plt

# Local modules
from fits import get_fits_header, get_fits_data
from fill_header import init_observatory 

def ron_gain_dark(my_instr="Mexman"):
    instrument = init_observatory(my_instr)
    gain = instrument['gain'] or 1
    ron = instrument['ron'] or 0
    dark_current = instrument['dark_current'] or 0
    log.info(f"Gain: {gain}, RON: {ron}, Dark current: {dark_current}")
        
    return ron, gain, dark_current

def detect_sources(image):
    ''' By Anna Marini
    Extract the light sources from the image
    '''
    # threshold = detect_threshold(image, nsigma=2.)
    # sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
    # kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    # kernel.normalize()

    if isinstance(image,str):
        image = get_fits_data(image)

    mask = make_source_mask(image, nsigma=2, npixels=5, dilate_size=11)
    mean,median,std = sigma_clipped_stats(image, sigma=3, mask=mask)
    daofind = DAOStarFinder(fwhm=3.0, threshold=5.*std)
    sources = daofind(image - median)
    # Pixel coordinates of the sources
    x = np.array(sources['xcentroid'])
    y = np.array(sources['ycentroid'])
    return x,y


def rescale(array):
    '''Take an array.  Rescale min to 0, max to 255, then change dtype,
    as opencv loves uint8 data type.  Returns the rescaled uint8 array.
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
                                    #frame='icrs',
                                    #unit="deg",
                                    radius=radius,
                                    catalog = 'Gaia', version=2)

    return catalog


def set_apertures(catalog, limit=16, r=10, r_in=15.5, r_out=25):
    '''From Anna Marini: get a catalog and
    set apertures and annulus for photometry.
    '''
    radec = catalog['ra','dec','phot_g_mean_mag']
    mask = radec['phot_g_mean_mag'] < limit
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
    
    ron,gain,dark_current = ron_gain_dark()

    for n in phot_table['residual_aperture_sum']:
        phot_table['S/N'] = final_sum / np.sqrt(final_sum
                                                + bkg_mean * pixar.area
                                                + ron
                                                + (gain/2)**2
                                                + dark_current)
    #phot_table['poisson_err'] = np.sqrt(final_sum)
    
    # Calculate errorbar
    # gain,ron = funzione che finisce con return gain,ron
    # if gain and ron:
    # signal_noise(pixar, pixan, gain, ron)
    # aggiungi colonna con l'errore.
    
    return(phot_table)


def apphot(filenames, reference=0, display=False, r=False, r_in=False, r_out=False):

    filenames = sorted(filenames)

    header0 = get_fits_header(filenames[reference])
    wcs0 = WCS(header0)
    
    catalog = load_catalog(wcs=wcs0)
    if r and r_in and r_out:
        apers = set_apertures(catalog, r=r, r_in=r_in, r_out=r_out)
    else:
        apers = set_apertures(catalog)

    tables = Table()

    if display:
        import pyds9
        d = pyds9.DS9("ds9")
        
    for filename in filenames:
        header = get_fits_header(filename)
        data = get_fits_data(filename)
        wcs = WCS(header)
    
        #catalog = load_catalog(wcs=wcs)
        #apers = set_apertures(catalog, r=r, r_in=r_in, r_out=r_out)
            
        phot_table = do_photometry(data, apers, wcs, obstime=header["MJD-OBS"])

        positions = SkyCoord(catalog['ra'], catalog['dec'],
                             frame='icrs',
                             unit=(u.deg,u.deg))                

        if display:
            d.set(f"file {filename}")
            
            # for i,pos in enumerate(positions):
            #     p = pos.to_pixel(wcs)
            #     circ = f'circle({p[0]}, {p[1]}, {10})'
            #     d.set("regions", circ)
            #     d.set("region", f"text {p[0]} {p[1]} "+"{"+str(i)+"}")
                
            for i,aper in enumerate(apers[0].to_pixel(wcs)) :
                circ = f'circle({aper.positions[0]}, {aper.positions[1]}, {aper.r})'
                d.set("regions", circ)
                d.set("region", f"text {aper.positions[0]}, {aper.positions[1]} "+"{"+str(i)+"}")
            
            for a in apers[1].to_pixel(wcs):
                circ = f'circle({aper.positions[0]}, {aper.positions[1]}, {aper.r_in})'
                d.set("regions", circ)
                circ = f'circle({aper.positions[0]}, {aper.positions[1]}, {aper.r_out})'
                d.set("regions", circ)
        
        tables.add_column(phot_table["residual_aperture_sum"], rename_duplicate=True)
            
        log.info(f"Done {filename}")

        return tables, phot_table

def plot(filenames):
    tables, phot_table = apphot(filenames, r=6, r_in=15.5, r_out=25)
    sources_number = len(tables)
    files_number =  len(tables[0])

    datas = [get_fits_header(f, fast=True)["MJD-OBS"] for f in filenames]
    tabellone = np.array([tables[k] for k in tables.keys() ])
    tabellone = np.insert(tabellone,  0, [datas], axis=1)
    t = tabellone[:,:1]
    flux = tabellone[:,1:files_number]
    ones = np.ones([files_number,sources_number])
    flux_err = phot_table['S/N']*ones
    
    fig1,ax = plt.subplots()
    n = list(range(sources_number))
    for number in n:
        ax.errorbar(t,flux[:,number],xerr= 0, yerr = flux_err[:,number],
                    fmt= ' ',
                    elinewidth = 2.5,
                    marker = 'o',markersize = 2)
    fig2,ax = plt.subplots()
    flux_1 = flux[:,1]
    err_1 = flux_err[:,1]
    ax.errorbar(t,flux_1,yerr=err_1,fmt= ' ', elinewidth = 2,
                ecolor = 'black', marker='o', markersize = 3.5) 
    
    flux_div = flux/flux[25] 
    err_div = np.sqrt(((flux_err/flux)**2))*flux_div
    
    fig3,ax = plt.subplots()
    for number in n:
        ax.errorbar(t,flux_div[:,number],yerr = err_div[:,number], fmt=' ',
                    elinewidth = 2, marker = 'o',
                    markersize = 2.5)

    #add check for saturated sources
    return(plt.show())
