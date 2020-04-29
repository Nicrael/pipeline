# !/usr/bin/env python3
# -*- coding: utf-8 -*-



####################################
# Graveyard of dead functions
####################################

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


# values1 = ['U', 'B', 'V']
# values2 = [1, 2]
# list(itertools.product(values1, values2))
# [('U', 1), ('U', 2), ('B', 1), ('B', 2), ('V', 1), ('V', 2)]

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
# datas = [ f.data for f in frames]
# ccds =  [ ccdp.CCDData(d, unit='adu') for d in datas ]
# ccds = [ ccdp.CCDData(get_fits_data(f.name), unit='adu') for f in frames ]

# path = Path('/home/dail/first/second/third')
# path.mkdir(parents=True, exist_ok=True)

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


#     pattern = glob.glob("gj3470/*/flat/*.fit*", recursive=True)
#     heads = [ r.get_fits_header(i) for i in pattern ]

#     key = 'FILTER'
#     values = { h[key] for h in heads } # distinct values

#     for value in values:
#         a = Time.now()

#         names = { p for p,h in zip(pattern, heads) if h[key] == value }
#         data = [r.get_fits_data(d) for d in names ]

#         data_norm = [ d/np.mean(d) for d in data ]
#         del data
#         dmaster = np.median(data_norm, axis=2)
#         del data_norm

#         print(f'{key} {value} -> {len(names)} elements.')
#         print(f'Done in {Time.now().unix - a.unix :.1f}s')

#     print(f'All done in {Time.now().unix - a.unix :.1f}s')
