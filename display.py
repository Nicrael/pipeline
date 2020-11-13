#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Display related functions.
'''

# System modules
import pyds9

# Local modules
import reduction as r


def ds9(*instance):
    '''
    Attach to a given ds9 instance or create a new one.
    '''
    if not instance:
        targets = pyds9.ds9_targets()
    else:
        targets = [str(instance[0])]
    dsn = pyds9.DS9(targets[0]) if targets else pyds9.DS9()
    return dsn


def show(*imgs, frame=1, target=False):
    '''
    Show or append a list of images or filenames in ds9.
    It is possible to choose a specific frame from which start to append.
    It is possible to Choose a specific ds9 target process.
    '''
    dsn = ds9() if target is False else ds9(target)
    dsn.set("tile yes")

    # # If some images are filenames, get their data first.
    lst = list(imgs)
    print(lst)
    for i, img in enumerate(lst):
        if isinstance(img, str):
            lst[i] = r.get_fits_data(img)
        imgs = tuple(lst)

    # If a list of fits is provided
    if str(frame) in ["first", "last", "prev", "next", "current"]:
        if frame == "current":
            frame = dsn.get("frame")
        dsn.set("frame "+frame)  # set to last
        frame = dsn.get("frame")  # get the id of last

    for i, img in enumerate(imgs, start=int(frame)):
        dsn.set("frame "+str(i))
        dsn.set_np2arr(img)


def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:]  # File(s). "1:" stands for "From 1 on".
    show(pattern)


# If called as a script
if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <list of FITS files>")
        sys.exit()

    main()
