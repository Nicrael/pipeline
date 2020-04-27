#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from multiprocessing import Process
# from pathlib import Path
# import dfits

#                  1               2            3                   4                 5           6
# biases ------> MBIAS
#        darks - MBIAS ->   darks_debiased
#                           darks_debiased -> MDARK
#        flats - MBIAS ->   flats_debiased
#                           flats_debiased  - MDARK -->  flats_debiased_dedarked
#                                                        flats_debiased_dedarked -> MFLAT
#      objects - MBIAS -> objects_debiased
#                         objects_debiased  - MDARK -> objects_debiased_dedarked
#                                                      objects_debiased_dedarked  / MFLAT -> objects_reduc
#

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

####################################
# shortcuts
####################################

import glob
import numpy as np
from astropy.time import Time

from reduction import *

#mask = oarpaf_mask(mbias, output_file=mout)
#reg = oarpaf_mask_reg(mask, output_file=f'{prod}-MASK-{keys}{value}.reg')

a = Time.now()

##########################################
# Data based
##########################################

# MBIAS
biases = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
mbias = combine(biases, method='median')

# MDARK
darks = glob.glob("gj3470/*/dark/*.fit*", recursive=True)
mdark = combine(darks, mbias=mbias, method='median')

# MFLAT
flats = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
mflat = combine(flats, mbias=mbias, mdark=mdark,
                normalize=True, method='median')

# CLEAN CUBE
obj_all = glob.glob("gj3470/*/object/*.fit*", recursive=True)
objects = dfits(obj_all).fitsort(['object']).unique_names_for(('GJ3470  ',))
clean_cube = combine(objects, mbias=mbias, mdark=mdark, mflat=mflat,
                     method='cube')

# CLEAN SLICES
for o in objects:
    clean_slice = combine(o, mbias=mbias, mdark=mdark, mflat=mflat, method='cube')

print(f'Done in {Time.now().unix - a.unix :.1f}s')



##########################################
# Filename+Data based
##########################################

# MBIAS
biases = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
generic(biases, keys=['ccdxbin'], method="median")

# MDARK
darks = glob.glob("gj3470/*/dark/*.fit*", recursive=True)
generic(darks, keys=['ccdxbin', 'exptime'], method="median",
                mbias='MBIAS.fits')

# MFLAT
flats = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
generic(flats, keys=['ccdxbin', 'filter'], method="median",
                mbias='MBIAS.fits', mdark='MDARK.fits', normalize=True)

# CLEAN CUBE
obj_all = glob.glob("gj3470/*/object/*.fit*", recursive=True)
objects = dfits(obj_all).fitsort(['object']).unique_names_for(('GJ3470  ',))

generic(objects, ['ccdxbin', 'filter'], method='cube',
                     mbias='MBIAS.fits', mdark='MDARK.fits', mflat="MFLAT.fits")

# CLEAN SLICES
generic(objects, ['ccdxbin', 'filter'], method='slice',
                     mbias='MBIAS.fits', mdark='MDARK.fits', mflat="MFLAT.fits")


####################################
# Main
####################################




def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".

    # dfits ~/desktop/oarpaf/test-sbig-stx/*.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'
    # dfits ~/desktop/oarpaf/test-sbig-stx/*.fits | fitsort IMAGETYP NAXIS1 DATE-OBS | grep 'Bias' |grep '4096' |grep '2019-11-22' | awk '{print $1}'

# 'abc' -> ['abc']



if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
