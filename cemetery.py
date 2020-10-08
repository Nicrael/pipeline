# !/usr/bin/env python3
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

####################################
# shortcuts
####################################


####################################
# Graveyard of dead functions
####################################

def mask(data, sigma=3, output_file=None, header=None):
    '''
    Create a bad pixel mask
    '''
    mask = sigma_clip(data, masked=True).mask.astype(int)
    if output_file:
        write_fits(mask, output_file, header=header)

    return mask


def mask_reg(data, sigma=3, output_file=None):
    '''
    Create a bad pixel region table
    '''

    y,x = np.where(data == True)
    p = np.repeat("point ", y.size)
    t = [p, x+1,y+1]

    table = Table(t, names=['# ','## ','###']) # bleah

    if output_file:
        ascii.write(table, output_file, overwrite=True)

    return table


def load(filenames):
    '''
    Return the datas of several filenames in a cube .
    '''
    sort(filenames)
    datas = np.array([ get_fits_data(f) for f in filenames ])

    # Collapsing array of cubes (3,30,100,100) -> (90,100,100)
    if len(datas.shape) > 3 :
        datas = datas.reshape(-1, *datas.shape[-2:])

    return datas.squeeze()

def is_keyval_in_header(head, key, val):
    '''
    Alternative for 'x' in header and header['x'] == 'y'
    '''
    return key in head and head[key] == val


def is_keyval_in_file(filename, key, val):
    '''
    Alternative for 'x' in fits and fits['x'] == 'y'
    '''
    header = get_fits_header(filename)
    return is_keyval_in_header(header, keyword, value)


def generic2(db, product=None, keys=None, mask=None,
             mbias=None, normalize=False):

    tab = db.table
    head = db.heads[0]
    if mask:
        the_mask = [mask in tt for tt in tab['FULLPATH'] ]
        tab = tab[the_mask]
    if keys:
        tab = tab.group_by(keys)

    for g in tab.groups: # in questo caso solo uno
        #print(g)
        files = list(g["FULLPATH"])
        master = combine(files, method='median', mbias=mbias,
                         normalize=normalize)
        #message = hist(__name__, f"Combining: {files}")
        #head["HISTORY"] = message
        text = dict(zip(keys, list(g[keys][0]))) if keys else None
        out = output_file(product=product, text=text)
        write_fits(data=master, output_file=out, header=head, fast=True)



import fitsio

def dfits(self, pattern, hdu=0):
    heads = [ fitsio.read_header(f, hdu) for f in pattern ]
    return self

def fitsort(self, keys, hdu=0):
    pattern = self.pattern
    heads = self.heads
    sort = [ (p,(tuple(h[k] for k in keys))) for p,h in zip(self.pattern,heads) ]
    self.keys = keys
    self.data = sort
    self.data = sort
    return self

def grep(self, value):
    sort = self.sort
    gr = { s[0] for s in sort if s[1] == value }
    self.data = gr
    return self

def unique(self, keys=None):
    if not keys:
        keys = self.keys
    heads = self.heads
    uns = { tuple(h[k] for k in keys) for h in heads }
    self.data = uns
    return self



def sethead2(head):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file
    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-test.json') as jfile:
        header_dict = json.load(jfile)# ['primary']

    template = fits.Header([tuple(d) for d in header_dict ] )

    diff = fits.HeaderDiff(template, head)

    for c in diff.common_keywords:
        template[c] = head[c]

    return template



def sethead3(head):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file
    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-dict.json') as jfile:
        header_dict = json.load(jfile)# ['primary']

    #card = [fits.Card(**c) for c in head_dict if key == c['keyword'] ][0]
    #form = '' if not card else card[0]["format"]

    for key in head :
        key =  key.lower()
        val = head[key]
        if key in header_dict: # keyword in the dictionary
            form = header_dict[key][0]  # Format in the dictionary
            comm =  header_dict[key][1] # Comment in the dictionary
            if 'd' in form:
                value = int(round(float(val)))
            elif 'b' in form:
                value = bool(val)
            elif 'f' in form:
                decimals = [ int(v) for v in form if v.isdigit() ]
                value = round(val, decimals[0])
            else: # 's' in form:
                value = f'{val:{form}}'
            print(key, val, "→", form, "→", value)
            head[key] = (value, comm)
        else:
            print(f'{key} not in dict')
            #value = val
            # comm = head.comments[key]

    return head




    def names(self, keys=None):
        if not keys:
            keys = self.keys
        heads = self.heads
        uns = { tuple(h[k] for k in keys) for h in heads }
        self.data = uns
        return self

    def unique(self, keys=None):
        if not keys:
            keys = self.keys
        heads = self.heads
        uns = { tuple(h[k] for k in keys) for h in heads }
        self.data = uns
        return self

    def grep(self, value):
        sort = self.data
        gr = { s[0] for s in sort if s[1] == value }
        self.data = gr
        return self


def asdasdasd():
    heads = [ get_fits_header(i, with_fitsio=with_fitsio) for i in pattern ]
    values = {tuple(h[k] for k in keys) for h in heads}
    # {('U', 10), ('U', 20), ('B', 10), ('B', 20)

    for value in values: # ('U', 10)

    # FILTER!!! WOW!
    names = { p for p,h in zip(pattern,heads) if tuple(h[k] for k in keys) == value }


def frame_dict(filename, with_data=False):
    '''
    Create a dictionary related to an observation frame
    '''
    fd = {
        'name' : filename,
        'head' : get_fits_header(filename),
        'data' : None,
        }
    if with_data:
        fd['data'] = get_fits_data(filename)

    return fd


class AttrDict(dict):
    '''
    Create an objects where properties are dict keys.
    '''
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def join_fits_header(pattern):
    '''
    Join the header of list of fits files in a list.
    '''
    heads = np.array([ get_fits_header(f) for f in pattern ])
    return heads


def join_fits_data(pattern):
    '''
    Join the data of a list of fits files in a tuple.
    Tuple format is useful for stacking in a data cube.
    '''
    datas = np.array([ get_fits_data(f) for f in pattern ])
    return datas


def stack_fits_data(datas):
    '''
    Stack a list of fits datas in a data cube.
    It is useful to perform pixel-per-pixel operations,
    such as an average.
    '''
    datacube = np.dstack(datas)
    return datacube


def median_datacube(datacube):
    '''
    Make a median of a data cube.
    '''
    datatype=datacube.dtype
    median = np.median(datacube, axis=2)
    return median.astype(datatype)


def average_datacube(datacube):
    '''
    Make an average of a data cube.
    '''
    datatype=datacube.dtype
    average = np.average(datacube, axis=2)
    return average.astype(datatype)


def oarpaf_combine(pattern, method='median', output_file=None, header=None):
    '''
    Custom master bias routine.
    Calculates the master bias of a list of of fits files.
    Default combining method is median.
    No output file is provided by default.
    '''
    joined_fits = join_fits_data(pattern)
    datacube = stack_fits_data(joined_fits)
    del joined_fits # saving memory
    if method is 'average':
        combined_data = average_datacube(datacube)
    else:
        combined_data = median_datacube(datacube)
    del datacube # saving memory

    header = get_fits_header(pattern[0]) if header else None

    if output_file:
        write_fits(combined_data, output_file, header=header)

    return combined_data

def new_header():
    return fits.PrimaryHDU().header

def to_list(arg):
    if type(arg) is not list: arg = [ arg ]
    return arg

def is_number(s):
    '''
    Check if a string contains a (float) number.
    Useful to test decimal or sexagesimal coordinates.
    '''
    try:
        float(s)
        return True
    except ValueError:
        return False


def to_number(s):
    try:
        return int(s)
    except ValueError:
        return float(s)


def get_fits_data_or_header(filename,get):
    '''
    Return the header or the data of a fits file.
    '''
    which_hdu = choose_hdu(filename)
    with fits.open(filename) as hdul:
        if get is 'header':
            return hdul[which_hdu].header;
        elif get is 'data':
            return hdul[which_hdu].data;
        else:
            return


def get_fits_data2(filename):
    '''
    Return the data of the fits file.
    Alternative method based on fitsio.
    '''
    which_hdu = choose_hdu(filename)
    with fitsio.FITS(filename) as f:
        data = f[which_hdu].read()
        return data


From nested dict (json) to object
class obj(object):
    def __init__(self, d):
        for a, b in d.items():
            if isinstance(b, (list, tuple)):
                setattr(self, a, [obj(x) if isinstance(x, dict) else x for x in b])
            else:
                setattr(self, a, obj(b) if isinstance(b, dict) else b)



def frame(filename):
    '''
    Create a frame object related to an observation frame
    '''
    fr = AttrDict(frame_dict(filename))

    return fr


def frame_list(pattern):
    '''
    Create a list of frames from filename pattern
    '''
    list1 = []
    for filename in pattern:
        list1.append(frame(filename))

    return list1


def combine_old(pattern, keys=[], method='average', normalize=False, fast=True, header=False):
    '''
    Combine a pattern of images using average (default) or median.
    Loops over a list of keywords. normalize=True to combine flats.
    '''

    print('Getting headers of all files in pattern')
    dlist = dfits(pattern).fitsort(keys)
    # [ (name1, ('U', 10)), (name2, ('U', 20)) ]

    for value in dlist.unique_values: # ('U', 10)
        a = Time.now()
        names = dlist.unique_names_for(value)

        datas = np.array([get_fits_data(d, fast=fast) for d in names ])
        datatype = datas.dtype

        if normalize:
            datas = np.array([ d/np.mean(d) for d in datas ])

        if method is 'average':
            combined = np.average(datas, axis=0).astype(datatype)
        else:
            combined = np.median(datas, axis=0).astype(datatype)
        del datas # saving memory

        print(f'{keys} {value} -> {len(names)} elements.')
        print(f'Done in {Time.now().unix - a.unix :.1f}s')

        write_fits(combined, 'MBIAS-test.fits')

        #return combined


def correct(pattern, master, keys=[], operation=None, fast=True, header=None):
    '''
    Take a pattern of file names. Correct data against a master.
    Use To subtract or divide.
    '''

    print('Getting headers of all files in pattern')
    dlist = dfits(pattern).fitsort(keys)
    # [ (name1, ('U', 10)), (name2, ('U', 20)) ]

    for value in dlist.unique_values: # ('U', 10)
        a = Time.now()
        names = dlist.unique_names_for(value)

        for name in names:
            if operation != 'flat':
                datas = get_fits_data(name, fast=fast) - master
            else:
                datas = get_fits_data(name, fast=fast) / master
            #yield

        print(f'{keys} {value} -> {len(names)} elements.')
        print(f'Done in {Time.now().unix - a.unix :.1f}s')


def subtract(pattern, master, keys=[], fast=True):
    '''
    Subtract two images
    '''
    if len(master) == 0: # Works for both [] and np.array
        master = 0
    return correct(pattern, master, keys=keys, fast=fast)


def divide(pattern, master, keys, fast=True):
    '''
    Divide two images
    '''
    if len(master) == 0: # Works for both [] and np.array
        master = 1
    return correct(pattern, master, keys=keys, fast=fast)


def ccdproc_mbias(pattern, method='median', output_file=None, header=None):

    joined_fits = [get_fits_data(p) for p in pattern]
    ccd_data_list = [ccdp.CCDData(d, unit="adu") for d in joined_fits ]

    combined_data = ccdp.combine(ccd_data_list, method=method, unit=u.adu,
                                 dtype=np.uint16, mem_limit=512e6, #)
                                 sigma_clip=True)  #, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 #sigma_clip_func=np.ma.median) #, signma_clip_dev_func=mad_std)
                                 #combined_bias.meta['combined'] = True

    if output_file:
        combined_bias.write(output_file, header=header, overwrite=True)

    return combined_data


values1 = ['U', 'B', 'V']
values2 = [1, 2]
list(itertools.product(values1, values2))
[('U', 1), ('U', 2), ('B', 1), ('B', 2), ('V', 1), ('V', 2)]

pattern
heads = [ get_fits_header(f) for f in pattern ]
sub_heads = is_keyval_in_header(heads, 'filter', 'vacio + B3')
frames = [ frame(f) for f in pattern if is_keyval_in_file(f, 'filter', 'vacio + B3') ]
datas = [ get_fits_data(f) for f in pattern ]
ccds =  [ ccdp.CCDData(d, unit='adu') for d in datas ]

pattern
all_frames = frame_list(pattern)
frames = [ f for f in all_frames if f.head['filter']  == 'vacio + B3']
files = [ f.name for f in frames]
heads = [ f.head for f in frames]
datas = [ f.data for f in frames]
ccds =  [ ccdp.CCDData(d, unit='adu') for d in datas ]
ccds = [ ccdp.CCDData(get_fits_data(f.name), unit='adu') for f in frames ]

path = Path('/home/dail/first/second/third')
path.mkdir(parents=True, exist_ok=True)

a = [1, 2, 3]
b = [4, 5, 6]
[list(zip(a, p)) for p in permutations(b)]
[[(1, 4), (2, 5), (3, 6)],
  [(1, 4), (2, 6), (3, 5)],
  [(1, 5), (2, 4), (3, 6)],
  [(1, 5), (2, 6), (3, 4)],
  [(1, 6), (2, 4), (3, 5)],
  [(1, 6), (2, 5), (3, 4)]]

a = ['U', 'B', 'V']
b = [1, 2]
list(itertools.product(a, b))
[('U', 1), ('U', 2), ('B', 1), ('B', 2), ('V', 1), ('V', 2)]


a = Time.now()
b = Time.now()
c = b.unix - a.unix


    pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
    heads = [ r.get_fits_header(i) for i in pattern ]

    key = 'FILTER'
    values = { h[key] for h in heads } # distinct values

    for value in values:
        a = Time.now()

        names = { p for p,h in zip(pattern, heads) if h[key] == value }
        data = [r.get_fits_data(d) for d in names ]

        data_norm = [ d/np.mean(d) for d in data ]
        del data
        dmaster = np.median(data_norm, axis=2)
        del data_norm

        print(f'{key} {value} -> {len(names)} elements.')
        print(f'Done in {Time.now().unix - a.unix :.1f}s')

    print(f'All done in {Time.now().unix - a.unix :.1f}s')

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

##########################################
# Generic2
##########################################

# INIT
pattern = "gj3470/*/*/*.fit*"
filenames = glob.glob(pattern, recursive=True)
db = minidb(filenames)


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
# With mini db table
####################################


db = table(filenames)
imagetyps = group(db, ["IMAGETYP", "FILTER"], "FULLPATH")

for s in  imagetyps:
    print(s, len(imagetyps[s]))


def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # File(s). "1:" stands for "From 1 on".


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :    # argv[0] is the filename.
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
