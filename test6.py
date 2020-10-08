# !/usr/bin/env python3
# -*- coding: utf-8 -*-


from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.table import vstack, Table
from astropy.coordinates import SkyCoord
from reduction import get_fits_header, get_fits_data
import glob
import numpy as np

from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats

pattern = "arp-data-2020-05-15T12:59:46/solved/*new"
filenames = sorted(glob.glob(pattern))

orig_datas = np.array([ get_fits_data(f, fast=True) for f in sorted(filenames) ])
wcss = np.array([ WCS(get_fits_header(f, fast=True)) for f in sorted(filenames) ])

#datas = orig_datas[:, 200:300, 530:630] 
#datas = orig_datas[:, 390:440, 570:620] 
datas = orig_datas[:, 300:380, 620:700] 

y,x=np.where(datas[0])

#c = SkyCoord.from_pixel(x,y, wcs=wcss[0])
#d = datas[0]

table = Table()

for i,d in enumerate(datas):
    c = SkyCoord.from_pixel(x,y, wcs=wcss[i])

    mask = make_source_mask(d, nsigma=2, npixels=5, dilate_size=11)
    mean,median,std = sigma_clipped_stats(d, sigma=3, mask=mask)
    
    v = d.reshape(d.size)
    rows = Table([c.ra, c.dec, v-median])
    table = vstack([rows,table])
    print(i) 

ascii.write(table, "rdv.dat") 

