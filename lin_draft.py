#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting

from fits import get_fits_header, get_fits_data

def blockshaped(arr, nrows=16, ncols=16):
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def counts(data):# dark_0.dat contiene i dati "grezzi" del dark a minore tempo di esposizione.
    blocks = blockshaped(data)
    tab_counts = Table()
    tab_counts['counts_avg'] = np.array([np.mean(n) for n in blocks.T]) # Transposed matrix (nota per me)
    return(tab_counts)


def time(filename):
    header = get_fits_header(filename)
    int_time = header['EXPTIME']
    tab_time = Table()
    tab_time['integration_time'] = [int_time]
    return(tab_time)


def linearity(filenames):
    filenames = sorted(filenames)
    counts_table = Table()
    time_table = Table()
    for filename in filenames:
        data = get_fits_data(filename)
        tab_counts = counts(data)
        tab_time = time(filename)
        counts_table.add_column(tab_counts['counts_avg'], rename_duplicate = True)
        time_table.add_column(tab_time['integration_time'], rename_duplicate = True)
    return(counts_table, time_table)


def plot(filenames, dark = True): #STX
    counts_table, time_table  = linearity(filenames)
    
    time = np.squeeze(np.array([time_table[k] for k in time_table.keys()]))
    counts = np.array([counts_table[k] for k in counts_table.keys()])

    if dark == True:
        fig1,ax1 = plt.subplots()
        fig2,ax2 = plt.subplots()        
        for n in counts.T:  #Transposed 
            sigma = np.std(n, dtype= np.float64)
            ax1.errorbar(time, n, yerr = sigma,
                        fmt = ' ',
                        marker = 'o',
                        markersize = 3,
                        elinewidth=1)
                       
            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('Counts')
            ax1.legend(['Block1', 'Block2', 'Block3', 'Block4', 'Block5', 'Block6',
                        'Block7', 'Block8', 'Block9', 'Block10', 'Block11',
                        'Block12', 'Block13', 'Block14', 'Block14', 'Block16'])

            norm = n/max(n)
            sigma_norm = np.std(norm, dtype=np.float64)
            ax2.errorbar(time, norm, yerr = sigma_norm,
                        fmt = ' ',
                        marker = 'o',
                        markersize = 3,
                        elinewidth=1)
               
            ax2.set_xlabel('Time (s)')
            ax2.set_ylabel('Counts')
            ax2.legend(['Block1', 'Block2', 'Block3', 'Block4', 'Block5', 'Block6',
                        'Block7', 'Block8', 'Block9', 'Block10', 'Block11',
                        'Block12', 'Block13', 'Block14', 'Block14', 'Block16'])

    else:
        plt.figure()
        fit = fitting.LinearLSQFitter()
        or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip, niter=3, sigma=2)
        line_init = models.Linear1D()
        
        for n in counts.T:
            fitted_line, mask = or_fit(line_init, time, n)
            filtered_data = np.ma.masked_array(n, mask=mask)
            plt.plot(time,n, 'ko', fillstyle='none')
            plt.plot(time, filtered_data, 'ko')
            plt.plot(time, fitted_line(time), 'k-')
            plt.xlabel('Time (s)')
            plt.ylabel('Counts')

    return(plt.show())


