
import glob
import numpy as np
from astropy.time import Time
from multiprocessing import Process
from pathlib import Path

import reduction as r


####################################
# Multi-process test
####################################

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
