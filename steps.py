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


# 111111111111111111111111111111111
# biases ------> MBIAS
# Loop over binning

prod = 'MBIAS'
keys =  ['CCDXBIN']
pattern = glob.glob("gj3470/*/bias/*.fit*", recursive=True)

r.combine(keys, pattern, prod, normalize=False)

# 2222222222222222222222222222222222
# darks - MBIAS ->   darks_debiased
# flats - MBIAS ->   flats_debiased
# objects - MBIAS -> objects_debiased
# Loop over binning

prod = 'debiased'
keys =  ['CCDXBIN']
master = 0

pattern = glob.glob("gj3470/*/d*/*.fit*", recursive=True)
r.subtract(keys, pattern, prod, master=None)
pattern = glob.glob("gj3470/*/f*/*.fit*", recursive=True)
r.subtract(keys, pattern, prod, master=None)
pattern = glob.glob("gj3470/*/o*/*.fit*", recursive=True)
r.subtract(keys, pattern, prod, master=None)

# pattern = glob.glob("gj3470/*/[dfo]*/*.fit*", recursive=True)
# r.subtract(keys, pattern, prod, master=None)

# 3333333333333333333333333333333333
# darks_debiased -> MDARK
# Loop over binning and exptime

prod = 'MDARK'
keys =  ['CCDXBIN', 'EXPTIME']
pattern = glob.glob("gj3470/*/dark/*.fit*", recursive=True)

r.combine(keys, pattern, prod, normalize=False)

# 4444444444444444444444444444444444
#   flats_debiased  - MDARK -->  flats_debiased_dedarked
# objects_debiased  - MDARK -> objects_debiased_dedarked
# Loop over binning and exptime

prod = 'dedarked'
keys =  ['CCDXBIN', 'EXPTIME']
pattern = glob.glob("gj3470/*/[fo]*/*.fit*", recursive=True)
master = 0

r.subtract(keys, pattern, prod, master=master)

# 5555555555555555555555555555555555
# flats_debiased_dedarked -> MFLAT
# Loop over binning and filter

prod = 'MFLAT'
keys =  ['CCDXBIN', 'FILTER']
pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)

r.combine(keys, pattern, prod, normalize=True)

# 6666666666666666666666666666666666
# objects_debiased_dedarked  / MFLAT -> objects_reduc
# Loop over binning and filter

prod = 'reduc'
keys =  ['CCDXBIN', 'FILTER']
pattern = glob.glob("gj3470/*/object/*.fit*", recursive=True)
master = 1

r.divide(keys, pattern, prod, master=master)


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
