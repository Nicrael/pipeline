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
       
    ron, gain, dark_current = ron_gain_dark()
    n_pix_aperture = pixan.area


    phot_table['S/N'] = final_sum / np.sqrt(final_sum
                                            + bkg_mean * pixar.area
                                            + ron
                                            + ((gain/2)**2)*n_pix_aperture
                                            + dark_current*n_pix_aperture)

##    phot_table['poisson_err'] = np.sqrt(final_sum)
                
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

def plot(filenames, flux_stable_ref = 1, limit = 65000, dat_file = 'gj3470-defot.dat'):
    filenames = sorted(filenames) 
    tables, phot_table = apphot(filenames, r=6, r_in=15.5, r_out=25)
    sources_number = len(tables)
    files_number =  len(tables[0])

    datas = [get_fits_header(f, fast=True)["MJD-OBS"] for f in filenames]
    tabellone = np.array([tables[k] for k in tables.keys()])
    tabellone = np.insert(tabellone,  0, [datas], axis=1)
    t = tabellone[:,:1]
    flux = tabellone[:,1:files_number]
    ones = np.ones([files_number,sources_number]) # in order to obtain consistent error matrix
    flux_err = phot_table['S/N']*ones
    magnitude = -2.5*np.log10(flux/flux_stable_ref) # with an ideal source who's flux == 1. 
    err_log = 0.434*(flux_err/flux) # for error given by logaritm base 10
    
    defot_table = ascii.read(dat_file)
##    defot_time_JD = Time(defot_table['col2'], format= 'jd')
##    defot_time = defot_time_JD.mjd # conversion from JD to MJD
    defot_time = defot_table['col2'] - 2400000 # why aren't the previous lines working?

    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()
    fig5,ax5 = plt.subplots()

    # Now let's plot defot and apphot magnitudes with different reference stars
    target_ref1 = (magnitude[:,2]-magnitude[:,0]) # Target - star number 2 
    target_mag1 = (defot_table['col8']-defot_table['col9'])
    error_defot1 = np.sqrt(defot_table['col14']**2 + defot_table['col15']**2)
    err_mag1 = np.sqrt(err_log[:,2]**2 + err_log[:,0]**2) # For magnitudes subtractions
    ax1.errorbar(defot_time, target_mag1, yerr = error_defot1, elinewidth = 1,
                marker = 'o', markersize = 1.5)
    ax1.errorbar(t,target_ref1, yerr = err_mag1,
                 elinewidth = 1, 
                marker = 'o', markersize = 1.5)
    ax1.legend(('Defot Plot', 'Apphot Plot'))
    ax1.set_xlabel('Time (MJD)')
    ax1.set_ylabel('Magnitude')
    ax1.set_title('Comparison using star 2')

    target_ref2 = (magnitude[:,2]-magnitude[:,1])
    target_mag2 = (defot_table['col8']-defot_table['col10'])
    error_defot2 = np.sqrt(defot_table['col14']**2 + defot_table['col16']**2)
    err_mag2 = np.sqrt(err_log[:,2]**2 + err_log[:,1]**2)
    ax2.errorbar(defot_time, target_mag2, yerr = error_defot2, elinewidth = 1,
                marker = 'o', markersize = 1.5)
    ax2.errorbar(t,target_ref2, yerr = err_mag2,
                 elinewidth = 1, 
                marker = 'o', markersize = 1.5)
    ax2.legend(('Defot Plot', 'Apphot Plot'))
    ax2.set_xlabel('Time (MJD)')
    ax2.set_ylabel('Magnitude')
    ax2.set_title('Comparison using star 3')
    

    target_ref3 = (magnitude[:,2]-magnitude[:,19])
    target_mag3 = (defot_table['col8']-defot_table['col11'])
    error_defot3 = np.sqrt(defot_table['col14']**2 + defot_table['col17']**2)
    err_mag3 = np.sqrt(err_log[:,2]**2 + err_log[:,19]**2)
    ax3.errorbar(defot_time, target_mag3, yerr = error_defot3, elinewidth = 1,
                marker = 'o', markersize = 1.5)
    ax3.errorbar(t,target_ref3, yerr = err_mag3,
                 elinewidth = 1, 
                marker = 'o', markersize = 1.5)
    ax3.legend(('Defot Plot', 'Apphot Plot'))
    ax3.set_xlabel('Time (MJD)')
    ax3.set_ylabel('Magnitude')
    ax3.set_title('Comparison using star 4')
    

    target_ref4 = (magnitude[:,2]-magnitude[:,14])
    target_mag4 = (defot_table['col8']-defot_table['col12'])
    error_defot4 = np.sqrt(defot_table['col14']**2 + defot_table['col18']**2)
    err_mag4 = np.sqrt(err_log[:,2]**2 + err_log[:,14]**2)
    ax4.errorbar(defot_time, target_mag4, yerr = error_defot4, elinewidth = 1,
                marker = 'o', markersize = 1.5)
    ax4.errorbar(t,target_ref4, yerr = err_mag4,
                 elinewidth = 1, 
                marker = 'o', markersize = 1.5)
    ax4.legend(('Defot Plot', 'Apphot Plot'))
    ax4.set_xlabel('Time (MJD)')
    ax4.set_ylabel('Magnitude')
    ax4.set_title('Comparison using star 5')
    

    target_ref5 = (magnitude[:,2]-magnitude[:,11])
    target_mag5 = (defot_table['col8']-defot_table['col13'])
    error_defot5 = np.sqrt(defot_table['col14']**2 + defot_table['col19']**2)
    err_mag5 = np.sqrt(err_log[:,2]**2 + err_log[:,11]**2)
    ax5.errorbar(defot_time, target_mag5, yerr = error_defot5, elinewidth = 1,
                marker = 'o', markersize = 1.5)
    ax5.errorbar(t,target_ref5, yerr = err_mag5,
                 elinewidth = 1, 
                marker = 'o', markersize = 1.5)
    ax5.legend(('Defot Plot', 'Apphot Plot'))
    ax5.set_xlabel('Time (MJD)')
    ax5.set_ylabel('Magnitude')
    ax5.set_title('Comparison using star 6')

    return(plt.show())
     
