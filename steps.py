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

from astropy import log
import glob
import numpy as np
from astropy.time import Time

from reduction import *
import naming
import fill_header
from sorters import dfits

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
objects = dfits(obj_all, fast=True).fitsort(['object']).unique_names_for(('GJ3470  ',))

clean_cube = combine(objects, mbias=mbias, mdark=mdark, mflat=mflat,
                     method='cube')

# CLEAN SLICES
for o in objects:
    clean_slice = combine(o, mbias=mbias, mdark=mdark, mflat=mflat)


##########################################
# Filename+Data based
##########################################

import glob
from reduction import master_bias, master_flat, correct_image
from naming import skeleton
from sorters import dfits

skeleton()

biases = glob.glob("gj3470/*/bias/*.fit*")
keys = ['ccdxbin']
master_bias(biases, keys)

flats = glob.glob("gj3470/*/flat/*.fit*")
keys = ['ccdxbin', 'filter']
mbias = 'arp.MBIAS.ccdxbin=2.fits'
master_flat(flats, keys, mbias=mbias)

obj_all = glob.glob("gj3470/*/object/*.fit*")
objects = dfits(obj_all, fast=True).fitsort(['object']).unique_names_for(('GJ3470  ',))
mflat = "arp.MFLAT.ccdxbin=2.filter=vacio+V3.fits"
correct_image(objects, keys, mbias=mbias, mflat=mflat)

solvituri = glob.glob("reduced/*new")
def solver(solvituri, header=False):
    import os

    cmd = f'solve-field \
    {pattern} \
    --crpix-center \
    --scale-units arcsecperpix \
    --scale-low 0.2 \
    --scale-high 0.6 \
    --ra 119.7433 \
    --dec 15.39145 \
    --radius 0.1 \
    --downsample 2 \
    --no-plots \
    --overwrite'

    os.system(cmd)


##########################################
# Generic2
##########################################

# INIT
pattern = glob.glob("gj3470/*/*/*.fit*", recursive=True)
db = minidb(pattern)


keys = ["CCDXBIN"]
mask = "bias"
product = "MBIAS"
combi(db, mask=mask, keys=keys, product=product)

keys = ["CCDXBIN", "FILTER"]
mask = "flat"
product = "MFLAT"
combi(db, mask=mask, keys=keys, product=product, normalize=True)




####################################
# All together
####################################

obj_all = glob.glob("gj3470/*/object/*.fit*", recursive=True)
objects = dfits(obj_all).fitsort(['object']).unique_names_for(('GJ3470  ',))

o = fill_header.init_observatory("Mexman")
o.filename = objects[0]
o.newhead()
new_header = o.nh
ooo = combine(objects[0], method='median',
              mbias="MBIAS.fits.fz", mflat="MFLAT-V.fits.fz"  )

name = naming.output_file(product="object")
write_fits(ooo, output_file=name, header=new_header)



####################################
# Main
####################################


from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, FK5, ICRS

# FITS
spm_ra = '8:00:26.29'     # penso il centro del campo, pixel 512
spm_dec = '+15:20:10.84'  # penso il centro del campo, pixel 516
# Scala in gradi: 0.000138 gradi per pixel.
# Il target si trova ai pixel 430, 385
# Differenza in pixel: -82,-131 → -0.0114,0.0182 gradi

spm_jd = 2458829.827152778
spm_equinox = 'J2019.9'

spm = SkyCoord(ra=spm_ra,
               dec=spm_dec,
               obstime=Time(spm_jd, format="jd"),
               frame="fk5",
               equinox=spm_equinox,
               unit=(u.hourangle, u.deg) )

spm_2000 = spm.transform_to(FK5(equinox="J2000"))

# SIMBAD
gj = SkyCoord.from_name("GJ3470").fk5

####################################
# With mini db table
####################################


db = table(pattern)
imagetyps = group(db, ["IMAGETYP", "FILTER"], "FULLPATH")

for s in  imagetyps:
    print(s, len(imagetyps[s]))




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
