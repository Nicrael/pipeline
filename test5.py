
import glob
import numpy as np
from astropy.time import Time
from multiprocessing import Process
from pathlib import Path

# dfits -x 1 gj3470/2019-12-1*/bias/* | fitsort zNAXIS1 zNAXIS2 JD


#path = Path('/home/dail/first/second/third')
#path.mkdir(parents=True, exist_ok=True)

pattern = glob.glob("gj3470/*/bias/*.fit*", recursive=True)

frames = [ r.frame(f) for f in pattern ]

key = 'CCDSUM'
values = [ f.head[key] for f in frames]

for unique_val in set(values):
    sub_frames = [ f for f in frames if f.head[key] == unique_val ]

    names = [ f.name for f in sub_frames ]

    # our method, average
    a = Time.now()
    mbias = r.oarpaf_combine(names, method="average")
    b = Time.now()
    c = b.unix - a.unix
    r.write_fits(data=mb, filename="oc.fits", header=frames[0].head)
    print(c)

    # # our method, median
    # a = Time.now()
    # mb = r.oarpaf_combine(names, method="median")
    # b = Time.now()
    # c = b.unix - a.unix
    # r.write_fits(data=mb, filename="ob.fits", header=frames[0].head)
    # print(c)

    # # ccdproc method, average, with clipping
    # a = Time.now()
    # mb = r.ccdproc_mbias(names, method="average")
    # b = Time.now()
    # c = b.unix - a.unix
    # r.write_fits(data=mb, filename="cc.fits", header=frames[0].head)
    # print(c)

    mask = r.oarpaf_mask(mb, output_file=f'MBIAS-{unique_val}.fits')
    r.oarpaf_mask_reg(mask, output_file='mask.reg')










pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)

heads = [ r.get_fits_header(f) for f in pattern ]

all_filters = [ h['filter'] for h in heads] #  if 'filter' in h]

unique_filters = set(all_filters)     # {'vacio + B3', 'vacio + I3', 'vacio + R3', 'vacio + U3', 'vacio + V3'}
# unique_filters = np.unique(filters) # array(['vacio + B3', 'vacio + I3', 'vacio + R3', 'vacio + U3', 'vacio + V3'], dtype='<U10')

pat = [ f for f in pattern if r.is_keyval_in_file(f, 'filter', 'vacio + I3') ]

a = Time.now()
mb = r.oarpaf_mbias(pat)
b = Time.now()
c = b.unix - a.unix

# a = Time.now()
# mb = r.ccdproc_mbias(pat)
# b = Time.now()
# c = b.unix - a.unix


def f(st, sleep):
    time.sleep(sleep)
    print(f"Ended   {st}")

def go2():
    for i in range(20):
        sleep = 0.2
        what = f"object_{i}"
        f(what, sleep)
        print(f"Started {what}")

def go3():
    for i in range(100):
        what = f"object_{i}"
        sleep = 0.1
        Process(target=f, args=(what, sleep) ).start()
        print(f"Started {what}")
