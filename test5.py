# !/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Process
import numpy as np


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


####################################
# Generators test
####################################


def many_to_one(inputs):
    output = np.average([d for d in inputs])
    return output

def many_to_many(inputs, master):

    for i in inputs:
        outputs = i - master
        yield outputs
