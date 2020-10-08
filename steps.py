#!/usr/bin/env python
# -*- coding: utf-8 -*-

#                  1               2            3                   4                 5           6
# biases ------> MBIAS
#        darks - MBIAS ->   darks_debiased -> MDARK
#        flats - MBIAS ->   flats_debiased  - MDARK --> flats_debiased_dedarked -> MFLAT
#      objects - MBIAS -> objects_debiased  - MDARK -> objects_debiased_dedarked / MFLAT -> objects_reduc
#

# biases, darks, flat
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract(*, MBIAS)
# 3 MDARK               = combine(darks_debiased)
# 4 *_debiased_dedarked = subtract(*_debiased, MDARK)
# 5 MFLAT               = combine(flats_debiased_dedarked)
# 6 objects_reduc       = divide(objects_debiased_dedarked, MFLAT)
#

# mask = mask(mbias, output_file=mout)
# reg = mask_reg(mask, output_file=f'{prod}-MASK-{keys}{value}.reg')

##########################################
# Data based
##########################################

# import glob
# from reduction import combine
# from sorters import dfits

# # MBIAS
# biases = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
# mbias = combine(biases, method='median')

# # MDARK
# darks = glob.glob("gj3470/*/dark/*.fit*", recursive=True)
# mdark = combine(darks, mbias=mbias, method='median')

# # MFLAT
# flats = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
# mflat = combine(flats, mbias=mbias, # mdark=mdark,
#                 normalize=True, method='median')

# # CLEAN CUBE
# obj_all = glob.glob("gj3470/*/object/*.fit*", recursive=True)
# objects = dfits(obj_all, fast=True).fitsort(['object']).unique_names_for(('GJ3470',))

# clean_cube = combine(objects, mbias=mbias, mflat=mflat, # mdark=mdark,
#                      method='cube')

# # CLEAN SLICES
# for o in objects:
#     clean_slice = combine(o, mbias=mbias, # mdark=mdark,
#                           mflat=mflat)


##########################################
# Filename+Data based
##########################################

import glob
from astropy.io import ascii
import numpy as np

from reduction import master_bias, master_flat, correct_image
from fits import get_fits_header
from naming import skeleton
from sorters import dfits
from fill_header import observatory, solver, init_observatory
from photometry import apphot
import glob

skeleton(date=True)

biases = glob.glob("gj3470/*/bias/*.fit*")
keys = ['ccdxbin']
master_bias(biases, keys)

flats = glob.glob("gj3470/*/flat/*.fit*")
keys = ['ccdxbin', 'filter']
mbias = 'arp.MBIAS.ccdxbin=2.fits'
master_flat(flats, keys, mbias=mbias)

obj_all = glob.glob("gj3470/*/object/*.fit*")
objects = dfits(obj_all, fast=True).fitsort(['object']).unique_names_for(('GJ3470',))
mflat = "arp.MFLAT.ccdxbin=2.filter=vacio+V3.fits"

correct_image(objects, keys, mbias=mbias, mflat=mflat, new_header="Mexman")

solver("*CLEAN*fits")

solved = sorted(glob.glob("arp-data-2020-05-19T20:06:35/solved/*CLEAN*.new")) 
datas = [get_fits_header(f, fast=True)["MJD-OBS"] for f in solved] 
airmass = [get_fits_header(f, fast=True)["AIRMASS"] for f in solved] 
tables = apphot(solved, r=6, r_in=15.5, r_out=25)
tabellone = np.array([tables[k] for k in tables.keys() ])
tabellone = np.insert(tabellone,  0, [datas, airmass], axis=1)
ascii.write( tabellone, "tabellone.txt", overwrite=True)



#plot f u ($1-58800):(-2.5*log10($5/($4+$14+$16))) w lp pt 7, g u ($2-2458800):($8+2.5*log10(10**(-$10*.4)+10**(-$12*.4)+10**(-$13*0.4) ))-0.000 w lp pt 7 lc rgb "orange"


####################################
# Main
####################################


# def main():
#     '''
#     Main function
#     '''
#     pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".

#     # dfits ~/desktop/oarpaf/test-sbig-stx/*.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'
#     # dfits ~/desktop/oarpaf/test-sbig-stx/*.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'

# # 'abc' -> ['abc']


# if __name__ == '__main__':
#     '''
#     If called as a script
#     '''
#     import sys

#     if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
#         print("Usage:  "+sys.argv[0]+" <list of FITS files>")
#         sys.exit()

#     main()
