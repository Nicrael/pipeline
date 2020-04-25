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

import reduction as r

#mask = r.oarpaf_mask(mbias, output_file=mout)
#reg = r.oarpaf_mask_reg(mask, output_file=f'{prod}-MASK-{keys}{value}.reg')

a = Time.now()

##########################################
# Data based
##########################################

# 1
biases = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
bbb = np.array([r.get_fits_data(b) for b in biases])
mbias = r.combine(bbb, method='median')

# 1,2,3
darks = glob.glob("gj3470/*/dark/*.fit*", recursive=True)
ddd = np.array([r.get_fits_data(d) for d in darks]).astype('uint16')
mdark = r.combine(ddd, mbias=mbias, method='median')

# 1,2,3,4,5
flats = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
fff = np.array([r.get_fits_data(f) for f in flats])
mflat = r.combine(fff, mbias=mbias, mdark=mdark, normalize=True, method='median')

# 1,2,3,4,5,6
obj_all = glob.glob("gj3470/*/object/*.fit*", recursive=True)
objects = r.dfits(obj_all).fitsort(['object']).unique_names_for(('GJ3470  ',))
for o in objects:
    ooo = r.get_fits_data(o)
    one = r.combine(ooo, mbias=mbias, mdark=mdark, mflat=mflat)

print(f'Done in {Time.now().unix - a.unix :.1f}s')


##########################################
# Filename+Data based
##########################################

# mbias
r.generic(biases, keys=['ccdxbin'])
# mdark
r.generic(darks, keys=['ccdxbin','exptime'], mbias=mbias)
# mflat
r.generic(flats, keys=['ccdxbin','filter'], mbias=mbias, mdark=mdark, normalize=True)


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


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
