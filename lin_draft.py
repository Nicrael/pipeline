#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from scipy.stats import linregress

from fits import get_fits_header, get_fits_data


def ins_temp(filename):
    header = get_fits_header(filename)
    instru = [header['INSTRUME']]
    if 'CCD-TEMP' in header:
        temp = [header['CCD-TEMP']]
    else:
        temp = [' ']

    return(instru, temp)


def counts(filename, bias1='1_CCD_Image_34.fits',
           bias2='1_CCD_Image_47.fits',
           bias3='1_CCD_Image_75.fits',
           bias4='2_CCD_Image_191.fits',
           bias5='2_CCD_Image_244.fits',
           bias6='2_CCD_Image_277.fits',
           bias7='3_CCD_Image_52.fits',
           bias8='3_CCD_Image_103.fits',
           bias9='3_CCD_Image_152.fits'):

    header = get_fits_header(filename)
    data = get_fits_data(filename)
    mean = np.mean(data)

    camera = header['INSTRUME']
    temp = header['CCD-TEMP']

    if camera == 'Atik Cameras':
        if temp<=0.3:
            mean_sub = mean-np.mean(get_fits_data(bias1))
        elif 0.3<temp<=1.3:
            mean_sub = mean-np.mean(get_fits_data(bias2))
        else:
            mean_sub = mean-np.mean(get_fits_data(bias3)) #intorno a 3
    elif camera == 'SBIG STL-11000 3 CCD Camera w/ AO':
        if -5.5<=temp<=-3.7:
            mean_sub = mean-np.mean(get_fits_data(bias4))
        elif temp<-5.5:
            mean_sub = mean-np.mean(get_fits_data(bias5))
        else:
            mean_sub = mean-np.mean(get_fits_data(bias6))
    elif camera == 'SBIG STX-16801 3 CCD Camera w/ AO':
        if temp<=-19.3:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias7))
        elif -19.2<temp<=-14.3:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias8))
        else:
            mean_sub = np.mean(get_fits_data(filename)[500:3500, 500:3500])-np.mean(get_fits_data(bias9))
            
    return(mean_sub)
    

def tab(filenames):
    filenames = sorted(filenames)
    time_table = Table()
    counts_unique_table = Table()
    for filename in filenames:
        header = get_fits_header(filename)
        counts_unique = counts(filename)
        time = header['EXPTIME']
        time_table.add_column([time], rename_duplicate = True)
        counts_unique_table.add_column([counts_unique], rename_duplicate = True)
    return(time_table, counts_unique_table)

def plot(filenames, Flat = True):
    time_table, counts_unique_table = tab(filenames)
    counts_unique = np.squeeze(np.array([counts_unique_table[k] for k in counts_unique_table.keys()]))
    time = np.squeeze(np.array([time_table[k] for k in time_table.keys()]))
    if Flat == True:      
        plt.figure()
        fit = fitting.LinearLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=2, sigma=2.5)
        line_init = models.Linear1D()

        fitted_line, mask = or_fit(line_init, time, counts_unique)
        filtered_data = np.ma.masked_array(counts_unique, mask=mask)
        plt.plot(time, counts_unique, 'ko', fillstyle='none')
        plt.plot(time, filtered_data, 'ko')
        plt.plot(time, fitted_line(time), 'k-')
        plt.xlabel('Time (s)')
        plt.ylabel('Counts')
        plt.title('Atik T=-1 Flat')
        plt.savefig('Atik_T=-1_Flat.pdf')

    else:    
        time = sorted(np.squeeze(np.array([time_table[k] for k in time_table.keys()])))

        slope, intercept, r_value, p_value, std_err = linregress(time, counts_unique)
        print(intercept)

        fig1,ax1 = plt.subplots()
        counts_norm = counts_unique/max(counts_unique)
        ax1.errorbar(time, counts_norm,
                     yerr = np.std(counts_norm, dtype=np.float64),
                    fmt = ' ',
                    marker = 'o',
                    markersize = 5,
                    elinewidth=1)
                   
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Counts')
        ax1.legend(['Block1', 'Block2', 'Block3', 'Block4', 'Block5', 'Block6',
                    'Block7', 'Block8', 'Block9', 'Block10', 'Block11',
                    'Block12', 'Block13', 'Block14', 'Block14', 'Block16'])
##            plt.title('Atik T=1 Dark')
##            plt.savefig('Atik_T=1_Dark.pdf')

##            sigma_norm = np.std(norm, dtype=np.float64)
           
    return(plt.show())
    
    
