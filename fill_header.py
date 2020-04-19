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
from reduction import get_fits_header, is_number, to_number, new_header


class observatory():

    def __init__(self, **kwargs):        
        '''
        Set default parameters.
        '''
                 
        # Header keywords
        self.ra = kwargs.get('ra', 'RA')
        self.dec = kwargs.get('dec', 'DEC')
        self.obj = kwargs.get('obj', 'OBJECT')
        self.obstime = kwargs.get('obstime', 'JD')
        self.binning = kwargs.get('binning')
        self.temperature = kwargs.get('temperature')
        self.pressure = kwargs.get('pressure')
        self.humidity = kwargs.get('humidity')
        
        # Constants        
        self.lon =  kwargs.get('lon')
        self.lat =  kwargs.get('lat')
        self.alt =  kwargs.get('alt', 0)
        self.unit =  kwargs.get('unit', (u.hour, u.deg))
        self.scale =  kwargs.get('scale')

        with open('cerbero-merged-array.json') as jfile:
            self.head_format = json.load(jfile)['primary']

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
        if self.ra and self.ra in head and self.dec in head:
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

        if self.binning and self.binning[0] in head:
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



    def sethead(self, h, key, val):
        ''' By Davide Ricci.
        Add a keyword in the header with comment and format
        taken from the json config file

        In [23]: "Hello %.2f" % 123.123
        Out[23]: 'Hello 123.12'

        In [30]: 'Hello, {:.2f}'.format(123.123)
        Out[30]: 'Hello, 123.12'
        '''
        
        card = [c for c in self.head_format if c["name"] == key]
            
        if not card:
            form = '{}'
            comm = None
        else:
            form = '{'+card[0]["format"]+'}'
            comm = card[0]["comment"]
        
        if is_number(val):
            value = form.format(val)
            value = float(value)
        else:
            value = val
  
        #print(key, val, "→", form, "→", value)
          
        if not card:
            h[key] = value
        else:
            h[key] = (value, card[0]["comment"])


    def newhead(self):
        nh = new_header()
        
        
        if not hasattr(self, 'plate'):
            self.detector()

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'w'):
            self.wcs()

        # location
        self.sethead(nh, "latitude", self.coord.location.lat.deg)
        self.sethead(nh, "longitud", self.coord.location.lon.deg)
        self.sethead(nh, "altitude", int(self.coord.location.height.to_value()) )

        # obstime
        self.sethead(nh, "mjd-obs", self.coord.obstime.mjd)
        self.sethead(nh, "jd", self.coord.obstime.jd)
        self.sethead(nh, "date-obs", self.coord.obstime.isot[:-4])
        self.sethead(nh, "st", self.coord.obstime.sidereal_time("mean").hour)
        self.sethead(nh, "equinox", self.coord.obstime.jyear_str)
        
        # detector
        self.sethead(nh, "plate", self.plate)

        # coord
        #self.sethead(nh, "ra", self.coord.ra.to_string(unit="hourangle",sep=":"))
        #self.sethead(nh, "dec", self.coord.dec.to_string(sep=":"))
        self.sethead(nh, "ra", self.coord.ra.deg)
        self.sethead(nh, "dec", self.coord.dec.deg)

        self.sethead(nh, "alt", self.coord.altaz.alt.deg)
        self.sethead(nh, "az", self.coord.altaz.az.deg)
        self.sethead(nh, "airmass", self.coord.altaz.secz.value)
        self.sethead(nh, "zdist", self.coord.altaz.zen.deg)
        
        sun_radec = get_sun(self.coord.obstime)
        sun_altaz = sun_radec.transform_to(self.coord.altaz)
        moon_radec = get_moon(self.coord.obstime)
        moon_altaz = moon_radec.transform_to(self.coord.altaz)

        self.sethead(nh, "sunalt", sun_altaz.alt.deg)
        self.sethead(nh, "sundist", sun_radec.separation(self.coord).deg)
        self.sethead(nh, "moonalt", moon_altaz.alt.deg)
        self.sethead(nh, "moondist", moon_radec.separation(self.coord).deg)

        self.sethead(nh, "filename", self.filename.split("/")[-1])
        
        # wcs
        #nh.extend(self.w.to_header(), update=True)
        
        #self.header.extend(w.to_header(), update=True)
        # self.header.rename_keyword("PC1_1", "CD1_1", force=True)

        nh.add_comment("Created by "+type(self).__name__)
        
        self.nh = nh

            
def main():    
    '''
    Main function
    '''

    with open('./instruments.json') as jfile:
        instruments = json.load(jfile)

    instrument = instruments[sys.argv[1]]
    o = observatory(**instrument)

    print(instrument)
    
    pattern = sys.argv[2:]
    for filename in pattern:
        o.filename = filename
        o.skycoord()
        o.newhead()
        print(o.nh['mjd-obs'])

        
if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 3 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <instrument as in json> <list of FITS files>")
        sys.exit()

    main()
