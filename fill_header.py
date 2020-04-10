#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.time import Time
import astropy.units as u
import numpy as np

# Our modules
from reduction import get_fits_header

class observatory():

    def __init__(self):
        '''
        Set default parameters.
        '''
        # No default filename or header
        self._filename = None
        self.header = None

    @property
    def filename(self):
        '''Laboratory image file name'''
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        self.header = get_fits_header(value)
        self._filename = value


    def config(self, filename=None):
        if filename is not None:
            self.filename()
            self.location()



    def location(self):

        #mexman
        longitud = self.header['LATITUDE']
        latitude = self.header['LONGITUD']
        altitude = self.header['ALTITUDE']

        #dfosc
        longitud = None
        latitude = None
        altitude = None

        #oarpaf
        longitud = None
        latitude = None
        altitude = None

        o = EarthLocation(lat=latitude*u.deg,
                          lon=longitud*u.deg,
                          height=altitude*u.m)


    def time(self):
        # all
        exptime = self.header['EXPTIME']

        #mexman
        obstime = self.header['JD'] # 2 sec difference wrt DATE...
        t = Time(obstime, format='jd')

        #dfosc
        obstime = self.header['MJD-OBS']
        t = Time(obstime, format='mjd')

        #oarpaf
        obstime = self.header['JD']
        t = Time(obstime, format='jd')

    def coordinates(self):
        #mexman
        ra  = self.header['RA']  # 18:56:10.8
        dec = self.header['DEC'] # 40:57:34.0
        ha  = self.header['HA']  # 02:54:11.3
        c = SkyCoord(ra=ra,
                     dec=dec,
                     unit=(u.hourangle, u.deg))

        #dfosc
        ra  = self.header['RA']  # 14.32572
        dec = self.header['DEC'] # -20.32528
        ha  = self.header['HA']  # 0.84511
        c = SkyCoord(ra=ra,
                     dec=dec,
                     unit=(u.deg, u.deg))

        #oarpaf
        ra  = None
        dec = None
        ha  = None


    def detector(self):

        #mexman
        binning = [self.header['CCDXBIN'],
                   self.header['CCDXBIN']]

        #mexman
        binning = list(map(int, self.header['CCDSUM'].split(' ')))

        #dfosc
        binning = [1, 1]

        #oarpaf
        binning = [self.header['XBINNING'],
                   self.header['YBINNING']]


    def meteo(self):

        #mexman
        temperature = self.header['XTEMP']*u.Celsius # 20*u.Celsius
        humidity = self.header['HUMIDITY']/100 # 0-1
        pressure = None    # 1000*u.hpa
        wavelength = None  # 550*u.nm

        #dfosc
        temperature = self.header['XTEMP']*u.Celsius # 20*u.Celsius
        humidity = self.header['HUMIDITY']/100 # 0-1
        pressure = None    # 1000*u.hpa
        wavelength = None  # 550*u.nm

        #oarpaf
        temperature = None
        humidity = None
        pressure = None    # 1000*u.hpa
        wavelength = None  # 550*u.nm


def fill_header(filename):

    observing_location = EarthLocation(lat=41.3*u.deg,
                                       lon=-74*u.deg,
                                       height=100*u.m)

    temperature = None # 20*u.Celsius
    pressure = None    # 1000*u.hpa
    humidity = None    # 0.2
    wavelength = None  # 550*u.nm

    observing_time = Time('2017-02-05 20:12:18')
    aa = AltAz(location=observing_location, obstime=observing_time)

    obj_radec = SkyCoord('4h42m', '-38d6m50.8s')
    obj_altaz = object_radec.transform_to(aa)
    zdist = object_altaz.zen # only in dfosc
    airmass = object_altaz.secz # only in mexman

    sun_radec = get_sun(observing_time)
    sun_altaz = sun_radec.transform_to(aa)

    moon_radec = get_moon(observing_time)
    moon_altaz = moon_radec.transform_to(aa)


def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # "1:" stands for "From 1 on".
    show(pattern)


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :
        print("Usage:  "+sys.argv[0]+" <Parameters>")
        sys.exit()

    main()
