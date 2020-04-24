import glob
import numpy as np
from astropy.time import Time

from multiprocessing import Process
from pathlib import Path

import reduction as r

import dfits

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

#mask = r.oarpaf_mask(mbias, output_file=mout)
#reg = r.oarpaf_mask_reg(mask, output_file=f'{prod}-MASK-{keys}{value}.reg')

biases = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
darks = glob.glob("gj3470/*/dark/*.fit*", recursive=True)
flats = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
objects = glob.glob("gj3470/*/object/*.fit*", recursive=True)

# 111111111111111111111111111111111
keys = ['CCDXBIN']
r.combine(biases, keys=keys)
mbias = []

# 2222222222222222222222222222222222
r.subtract(darks, mbias, keys=keys)
r.subtract(flats, mbias, keys=keys)
r.subtract(objects, mbias, keys=keys)
darks_debiased = []
flats_debiased = []
objects_debiased = []

# 3333333333333333333333333333333333
keys = ['CCDXBIN','EXPTIME']
r.combine(darks_debiased, keys=keys)
mdark = []

# 4444444444444444444444444444444444
keys = ['CCDXBIN','EXPTIME']
r.subtract(flats_debiased, mdark, keys=keys)
r.subtract(objects_debiased, mdark, keys=keys)
flats_debiased_dedarked = []
objects_debiased_dedarked = []

# 5555555555555555555555555555555555
keys = ['CCDXBIN','EXPTIME']
r.combine(flats_debiased_dedarked, keys=keys, normalize=True)
mflat = []

# 6666666666666666666666666666666666
keys = ['CCDXBIN', 'FILTER']
r.divide(objects_debiased_dedarked, mflat, keys=keys)

####################################

def single_master(pattern, keys=[], normalize=False):
    r.combine(pattern, keys=keys, normalize=normalize)


def single_reduc(filename, mbias=mbias, mdark=mdark, mflat=mflat, keys=[]):
    r.subtract(keys, filename, master=mbias)
    r.subtract(keys, filename_debiased, master=mdark)
    r.divide(keys, filename_debiased_dedarked, master=mflat)


def multiple_reduc(pattern, mbias, mdark, mflat):
    r.subtract(['CCDXBIN'], pattern, master=mbias)
    r.subtract(['CCDXBIN', 'EXPTIME'], pattern_debiased, master=mdark)
    r.divide(['CCDXBIN', 'FILTER'], pattern_debiased_dedarked, master=mflat)



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
