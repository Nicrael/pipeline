#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import json

# from astropy.io import fits

# Our modules
from reduction import get_fits_header, is_number, to_number


class headinfo():

    def __init__(self, ra="RA", dec='DEC', obj='OBJECT',
                 unit=(u.hour, u.deg), obstime='JD',
                 lon=None, lat=None, alt=0,
                 binning=None, scale=None,
                 temperature=None, pressure=None, humidity=None):        
        '''
        Set default parameters.
        '''
        
        # Header keywords
        self.ra = ra
        self.dec = dec
        self.obj = obj
        self.obstime = obstime
        self.binning = binning
        self.temperature = temperature
        self.pressure = pressure
        self.humidity = humidity
        
        # Constants        
        self.lon = lon
        self.lat = lat
        self.alt = alt
        self.unit = unit
        self.scale = scale

    def test(self, filename):
        '''
        Runs all methods to have all parameters.
        Mainly for test purposes
        '''
        self.filename = filename
            
        self.skycoord()
        self.detector()
        self.wcs()

    @property
    def filename(self):
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        self._filename = value
        self.head = get_fits_header(value)
        

    def skycoord(self):
        head = self.head

        # earthlocation
        location = EarthLocation(lon=self.lon,
                                 lat=self.lat,
                                 height=self.alt)
        
        # time
        timekey = head[self.obstime]        
        if 'MJD' in self.obstime:
            obstime = Time(timekey, format='mjd')
        elif self.obstime == 'JD':
            obstime = Time(timekey, format='jd')
        # elif isinstance(self.obstime, list):
        #     obstime = Time(Time(timekey[0]).unix+timekey[1])
        elif 'DATE' in self.obstime:
            obstime = Time(timekey)
        else:
            obstime = Time(timekey)

        obstime.location=location
        
        # skycoord
        if self.ra in head and self.dec in head:
            coord = SkyCoord(ra=head[self.ra],
                             dec=head[self.dec],
                             unit=self.unit)
        elif self.obj in head:
            coord = SkyCoord.from_name(head[self.obj])
        else:
            return

        # Meteo 
        if self.temperature and self.temperature in head:
            coord.temperature = head[self.temperature] *u.deg_C
        if self.pressure and self.pressure in head:
            coord.pressure = head[self.pressure] *u.hPa
        if self.humidity and self.humidity in head:
            coord.relative_humidity = head[self.humidity]/ 100 
        
        coord.location=location
        coord.obstime=obstime

        self.coord = coord
        return coord

    
    def detector(self):
        '''
        Manage keywords related to the detector.
        '''
        head = self.head

        if self.binning and self.binning in head:
            bins = [self.head[self.binning[0]],
                    self.head[self.binning[1]]]
        else:
            bins = [1, 1]

        if not self.scale:
            self.scale = 1

        plate = self.scale * bins[0]
        
        self.bins = bins
        self.plate = plate
        return plate

        
    def wcs(self):
        '''From Anna Marini.
        Provide WCS keywords to convert pixel coordinates of the
        files to sky coordinates. It uses the rotational matrix
        obtained in previous function (which_instrument)
        '''
        head = self.head

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'plate'):
            self.detector()

        angle = 0
        flip = 1
        if 'INSTRUME' in head and head['INSTRUME'] == 'Mexman':
            angle = np.pi/2
            flip = -1

        cd = self.plate*np.array([[np.cos(angle), np.sin(angle)*flip],
                                  [np.sin(angle), np.cos(angle)]])

        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = cd
        w.wcs.crval = [self.coord.ra.deg,
                       self.coord.dec.deg]
        w.wcs.crpix = [head['NAXIS1']/2,
                       head['NAXIS2']/2]

        # For other methods
        self.w = w
        return w

    
mexman = headinfo(lon=-115.4869,
                  lat=31.0292,
                  alt=2790)

dfosc = headinfo(lon="-70:44:14.7926",
                 lat="-29:15:27",
                 obstime="MJD-OBS")

sbig = headinfo(lon=9.2034,
                lat=44.5911,
                alt=1480)
