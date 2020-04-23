
import glob
import numpy as np
from astropy.time import Time

from multiprocessing import Process
from pathlib import Path
import itertools

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

r.divide(keys, pattern, prod)



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




    # values1 = ['U', 'B', 'V']
    # values2 = [1, 2]
    # list(itertools.product(values1, values2))
    # [('U', 1), ('U', 2), ('B', 1), ('B', 2), ('V', 1), ('V', 2)]





if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()



    #     ####################################
    #     cimitero
    #     ####################################

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

# a = [1, 2, 3]
# b = [4, 5, 6]
# [list(zip(a, p)) for p in permutations(b)]
# [[(1, 4), (2, 5), (3, 6)],
#   [(1, 4), (2, 6), (3, 5)],
#   [(1, 5), (2, 4), (3, 6)],
#   [(1, 5), (2, 6), (3, 4)],
#   [(1, 6), (2, 4), (3, 5)],
#   [(1, 6), (2, 5), (3, 4)]]

# a = ['U', 'B', 'V']
# b = [1, 2]
# list(itertools.product(a, b))
# [('U', 1), ('U', 2), ('B', 1), ('B', 2), ('V', 1), ('V', 2)]


# a = Time.now()
# b = Time.now()
# c = b.unix - a.unix


    # pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
    # heads = [ r.get_fits_header(i) for i in pattern ]

    # key = 'FILTER'
    # values = { h[key] for h in heads } # distinct values

    # for value in values:
    #     a = Time.now()

    #     names = { p for p,h in zip(pattern, heads) if h[key] == value }
    #     data = [r.get_fits_data(d) for d in names ]

    #     data_norm = [ d/np.mean(d) for d in data ]
    #     del data
    #     dmaster = np.median(data_norm, axis=2)
    #     del data_norm

    #     print(f'{key} {value} -> {len(names)} elements.')
    #     print(f'Done in {Time.now().unix - a.unix :.1f}s')

    # print(f'All done in {Time.now().unix - a.unix :.1f}s')
