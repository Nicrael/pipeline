
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
#                           flats_debiased  - MDARK -->  flats_debiased_dedarked -> MFLAT
#      objects - MBIAS -> objects_debiased
#                         objects_debiased  - MDARK -> objects_debiased_dedarked  / MFLAT -> objects_reduc
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

#path = Path('/home/dail/first/second/third')
#path.mkdir(parents=True, exist_ok=True)

####################################
# 1111111111111111111111111111111111
####################################

a = Time.now()

pattern = glob.glob("gj3470/*/bias/*.fit*", recursive=True)

prod = 'BIAS'

frames = [ r.frame(f) for f in pattern ]
key = 'CCDXBIN'
values = [ f.head[key] for f in frames ]

for unique_val in set(values):
    sub_frames = [ f for f in frames if f.head[key] == unique_val ]

    names = [ f.name for f in sub_frames ]

    out = f'MBIAS-{key}{unique_val}.fits'
    ref_header = frames[0].head

    a = Time.now()

    mbias = r.oarpaf_combine(names, method="average", output_file=out, header=ref_header)
    #mbias = r.oarpaf_combine(names, method="median", output_file=out, header=ref_header)
    #mbias = r.ccdproc_combine(names, method="average", output_file=out, header=ref_header)

    mout = f'{prod}-{key}{unique_val}.fits'
    mask = r.oarpaf_mask(mb, output_file=mout)
    reg = r.oarpaf_mask_reg(mask, output_file=f'{prod}-{key}{unique_val}-MASK.reg')

    b = Time.now()
    c = b.unix - a.unix
    print(c)


####################################
# 2222222222222222222222222222222222
####################################

a = Time.now()

pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)

frames = [ r.frame(f) for f in pattern ]
key = 'CCDXBIN'
values = [ f.head[key] for f in frames ]

def deb(f):
    raw_name = f.name
    raw_data = r.get_fits_data(f.name)
    ccd_debiased = raw_data - master_data

    b = Time.now()
    c = b.unix - a.unix

    return [raw_name, master_name, c]

def go4():
    for unique_val in set(values):
        sub_frames = [ f for f in frames if f.head[key] == unique_val ]

        master_name = f'MBIAS-{key}{unique_val}.fits'
        master_data = r.get_fits_data(master_name)

        for f in sub_frames:
            sss=deb(f)
            #sss=Process(target=deb(f) ).start()
            print(sss)

####################################
# 3333333333333333333333333333333333
####################################

####################################
# 4444444444444444444444444444444444
####################################

####################################
# 5555555555555555555555555555555555
####################################

####################################
# 6666666666666666666666666666666666
####################################

def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".

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
