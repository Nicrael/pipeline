#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.io import ascii

from fits import get_fits_header, get_fits_data

def instrument_temperature_time_dateobs(filename):
    header = get_fits_header(filename)
    instru = [header['INSTRUME']]
    if 'CCD-TEMP' in header:
        temp = [header['CCD-TEMP']]
    else:
        temp = [' ']
    time = [header['EXPTIME']]
    date_obs = header['DATE-OBS']
    return(instru, temp, time, date_obs)

def counts(filename,
           bias1='1_CCD_Image_34.fits',
           bias2='1_CCD_Image_47.fits',
           bias3='1_CCD_Image_75.fits',
           bias4='2_CCD_Image_191.fits',
           bias5='2_CCD_Image_244.fits',
           bias6='2_CCD_Image_277.fits',
           bias7='3_CCD_Image_52.fits',
           bias8='3_CCD_Image_103.fits',
           bias9='3_CCD_Image_152.fits'):
    
    header = get_fits_header(filename)

    camera = header['INSTRUME']
    temp = header['CCD-TEMP']

    if camera == 'Atik Cameras':
        if temp<=0.3:
            mean_sub = np.mean(get_fits_data(filename)[1007:3007, 631:2040])-np.mean(get_fits_data(bias1))
        elif 0.3<temp<=1.3:
            mean_sub = np.mean(get_fits_data(filename)[1007:3007, 631:2040])-np.mean(get_fits_data(bias2))
        else:
            mean_sub = np.mean(get_fits_data(filename)[1007:3007, 631:2040])-np.mean(get_fits_data(bias3)) #intorno a 3
    elif camera == 'SBIG STL-11000 3 CCD Camera w/ AO':
        if -5.5<=temp<=-3.7:
            mean_sub = np.mean(get_fits_data(filename)[500:3508, 500:2172])-np.mean(get_fits_data(bias4))
        elif temp<-5.5:
            mean_sub = np.mean(get_fits_data(filename)[500:3508, 500:2172])-np.mean(get_fits_data(bias5))
        else:
            mean_sub = np.mean(get_fits_data(filename)[500:3508, 500:2172])-np.mean(get_fits_data(bias6))
    elif camera == 'SBIG STX-16801 3 CCD Camera w/ AO':
        if temp<=-19.3:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias7))
        elif -19.2<temp<=-14.3:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias8))
        else:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias9))
            
    return(mean_sub)


def tables(filenames):
    tab_time = Table()
    tab_counts = Table()
    tab_temp = Table()
    tab_cam = Table()
    tab_date_obs = Table()
    tab_filename = Table()
    for filename in filenames:
        mean = [counts(filename)]
        camera, temp, exp_time, date_obs = instrument_temperature_time_dateobs(filename)
        t = Time(date_obs)
        jd = t.jd
        
        tab_time.add_column([exp_time], rename_duplicate = True)
        tab_counts.add_column([mean], rename_duplicate = True)
        tab_temp.add_column([temp], rename_duplicate = True)
        tab_cam.add_column([camera], rename_duplicate = True)
        tab_date_obs.add_column([jd], rename_duplicate = True)
        tab_filename.add_column([filename], rename_duplicate = True)
    
    return(tab_time, tab_counts, tab_temp, tab_cam, tab_date_obs, tab_filename)

def tab_ascii(filenames):
    filenames = sorted(filenames)
    time, counts, temp, cam, jd, filename = tables(filenames)
    time_arr = np.squeeze(np.array([time[k] for k in time.keys()]))
    counts_arr = np.squeeze(np.array([counts[k] for k in counts.keys()]))
    temp_arr = np.squeeze(np.array([temp[k] for k in temp.keys()]))
    cam_arr = np.squeeze(np.array([cam[k] for k in cam.keys()]))
    jd_arr = np.squeeze(np.array([jd[k] for k in jd.keys()]))
    file_arr = np.squeeze(np.array([filename[k] for k in filename.keys()]))


    tabascii = Table([file_arr, time_arr, jd_arr, counts_arr, temp_arr, cam_arr],
                     names=['file name', 'time (s)', 'time (jd)', 'counts (mean)',
                            'temperature (°C)', 'camera'])
    ascii.write(tabascii, 'Flat_sub_values',
                format='fixed_width',
                overwrite=True)
    return(tabascii)
