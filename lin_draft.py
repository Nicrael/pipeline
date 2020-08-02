#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt

from fits import get_fits_header, get_fits_data

def blockshaped(arr, nrows=16, ncols=16):
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def counts(data, mdark = 'dark_0.dat'): # dark_0.dat contiene i dati "grezzi" del dark a minore tempo di esposizione.
    mdark = ascii.read(mdark)
    mdark_data = np.array([mdark[k] for k in mdark.keys()])
    
    blocks_mdark = blockshaped(mdark_data)
    blocks = blockshaped(data)
    difference = blocks - blocks_mdark # sottraggo a ogni blocco il suo corrispettivo blocco dark_0

    tab_counts = Table()
    tab_counts['counts_avg'] = np.array([np.mean(n) for n in difference.T]) # Transposed matrix (nota per me)
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


def plot(filenames): #STX
    counts_table, time_table  = linearity(filenames)
    
    time = np.squeeze(np.array([time_table[k] for k in time_table.keys()]))
    counts = np.array([counts_table[k] for k in counts_table.keys()])

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
        ax1.set_ylabel('Counts (sub Dark1.fits)')
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
        ax2.set_ylabel('Normalized Counts (sub Dark1.fits)')
        ax2.legend(['Block1', 'Block2', 'Block3', 'Block4', 'Block5', 'Block6',
                    'Block7', 'Block8', 'Block9', 'Block10', 'Block11',
                    'Block12', 'Block13', 'Block14', 'Block14', 'Block16'])
# Ho assunto che i blocchi siano in ordine, ovvero da sinistra a destra e dall'alto in basso.
# A favore di questa ipotesi, ho visto che il primo elemento di "data" e il primo elemento di blocks[0]
# coincidono. Se ciò fosse vero, l'elemento più in alto sarebbe di un quadrato centrale, il 10.
# Mi pare che abbia senso, dovrebbero essere i quadrati dal flusso maggiore, poiché negli angolini si ha
# sempre un flusso minore.
# Ho normalizzato in vista del calcolo della DC, ma con i dati negativi chiaramente l'intercetta
# sarà negativa, quindi ho qualche dubbio sulla validità del risultato: ha senso, con questo andamento, dare
# una stima della corrente? 

    return(plt.show())


