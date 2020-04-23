
import glob
import numpy as np
from astropy.time import Time
from multiprocessing import Process
from pathlib import Path

import reduction as r

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
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 3 MDARK               = combine(darks_debiased)
# 4 *_debiased_dedarked = subtract_dark(*_debiased, MDARK)
# 5 MFLAT               = combine(flats_debiased_dedarked)
# 6 objects_reduc       = flat_correct(objects_debiased_dedarked, MFLAT)
#

# biases,      , flat
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 5 MFLAT               = combine(flats_debiased)
# 6 objects_reduc       = flat_correct(objects_debiased, MFLAT)
#

# biases, darks
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
# 3 MDARK               = combine(darks_debiased)
# 4 *_debiased_dedarked = subtract_dark(*_debiased, MDARK)
#

#         darks, flat
#
# 3 MDARK               = combine(darks)
# 4 *_dedarked          = subtract_dark(*, MDARK)
# 5 MFLAT               = combine(flats_dedarked)
# 6 objects_reduc       = flat_correct(objects_dedarked, MFLAT)
#

# biases
#
# 1 MBIAS               = combine(biases)
# 2 *_debiased          = subtract_bias(*, MBIAS)
#

#         darks
#
# 3 MDARK               = combine(darks)
# 4 *_dedarked          = subtract_dark(*, MDARK)
#

#                 flat
#
# 5 MFLAT               = combine(flats)
# 6 objects_reduc       = flat_correct(objects, MFLAT)
#

####################################
# shortcuts
####################################

# pattern
# heads = [ get_fits_header(f) for f in pattern ]
# sub_heads = is_keyval_in_header(heads, 'filter', 'vacio + B3')
# frames = [ frame(f) for f in pattern if is_keyval_in_file(f, 'filter', 'vacio + B3') ]
# datas = [ get_fits_data(f) for f in pattern ]
# ccds =  [ ccdp.CCDData(d, unit='adu') for d in datas ]

# pattern
# all_frames = frame_list(pattern)
# frames = [ f for f in all_frames if f.head['filter']  == 'vacio + B3']
# files = [ f.name for f in frames]
# heads = [ f.head for f in frames]
#### datas = [ f.data for f in frames]
# ccds =  [ ccdp.CCDData(d, unit='adu') for d in datas ]
# ccds = [ ccdp.CCDData(get_fits_data(f.name), unit='adu') for f in frames ]

#path = Path('/home/dail/first/second/third')
#path.mkdir(parents=True, exist_ok=True)

# a = Time.now()
# b = Time.now()
# c = b.unix - a.unix


####################################
# 1111111111111111111111111111111111
####################################

# biases ------> MBIAS
# Loop over binning

def loop1():

    prod = 'MBIAS'
    pattern = glob.glob("gj3470/*/bias/*.fit*", recursive=True)
    frames = [ r.frame(f) for f in pattern ]
    key = 'CCDXBIN'
    values = [ f.head[key] for f in frames ]

    for unique_val in set(values):

        sub_frames = [ f for f in frames if f.head[key] == unique_val ]
        names = [ f.name for f in sub_frames ]
        ref_header = frames[0].head

        out = f'{prod}-{key}{unique_val}.fits'

        mbias = r.oarpaf_combine(names, method="average", output_file=out, header=ref_header)
        #mbias = r.oarpaf_combine(names, method="median", output_file=out, header=ref_header)
        #mbias = r.ccdproc_combine(names, method="average", output_file=out, header=ref_header)

        mout = f'{prod}-MASK-{key}{unique_val}.fits'
        mask = r.oarpaf_mask(mb, output_file=mout)
        reg = r.oarpaf_mask_reg(mask, output_file=f'{prod}-MASK-{key}{unique_val}.reg')

        print(c)

####################################
# 2222222222222222222222222222222222
####################################

# darks - MBIAS ->   darks_debiased
# flats - MBIAS ->   flats_debiased
# objects - MBIAS -> objects_debiased

# Loop over binning

def loop2():

    prod = 'DEBIASED'
    pattern = glob.glob("gj3470/*/object/*.fit*", recursive=True)
    frames = [ r.frame(f) for f in pattern ]
    key = 'CCDXBIN'
    values = [ f.head[key] for f in frames ]

    for unique_val in set(values):
        sub_frames = [ f for f in frames if f.head[key] == unique_val ]

        out = f'{prod}-{key}{unique_val}.fits'
        master_name = f'{MBIAS}-{key}{unique_val}.fits'
        master_data = r.get_fits_data(master_name)

        for f in sub_frames:
            sss = step2(f)
            #sss = Process(target=deb(f) ).start()
            print(sss)

def step2(f):

    raw_name = f.name
    raw_data = r.get_fits_data(f.name)
    ccd_debiased = raw_data - master_data

    return [raw_name, master_name, c]

####################################
# 3333333333333333333333333333333333
####################################

# darks_debiased -> MDARK
# Loop over binning and exptime

def loop3():
    pass

####################################
# 4444444444444444444444444444444444
####################################

#   flats_debiased  - MDARK -->  flats_debiased_dedarked
# objects_debiased  - MDARK -> objects_debiased_dedarked

# Loop over binning and exptime

def loop4():
    pass

####################################
# 5555555555555555555555555555555555
####################################

# flats_debiased_dedarked -> MFLAT

# Loop over binning and filter

def loop5():
    pass

####################################
# 6666666666666666666666666666666666
####################################

# objects_debiased_dedarked  / MFLAT -> objects_reduc

def loop6():
    pass

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
