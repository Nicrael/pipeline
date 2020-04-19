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

# Our modules
from reduction import get_fits_header, is_number, to_number


class observatory():

    def __init__(self, filename=None):
        '''
        Set default parameters.
        '''
        # No default filename, header or instrument
        self._filename = None
        self.original = None # original header
        self.header = None
        self.params = None

        with open('instruments.json') as jfile:
            self.instruments = json.load(jfile)

        with open('cerbero-merged-array.json') as jfile:
            self.head_format = json.load(jfile)['primary']

        if filename is not None:
            self.filename = filename


    def config(self, filename=None):
        '''
        Runs all methods to have all parameters.
        Mainly for test purposes
        '''
        if filename is not None:
            self.filename = filename
            
        self.earthlocation()
        self.time()
        self.skycoord()
        self.detector()
        self.meteo()
        self.wcs()
        self.newhead()

    @property
    def filename(self):
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        self._filename = value
        self.original =  get_fits_header(value)
        self.header = self.original.copy()
        self.params = self.instruments[self.header['INSTRUME']]


    def earthlocation(self):
        '''By Davide Ricci.
        Manage keywords related to local parameters.
        '''

        param = self.params['location']
                
        location = EarthLocation(lon=param[0],
                                 lat=param[1],
                                 height=int(param[2]))

        # For other methods
        self.location = location
        return location

    def time(self):
        '''By Anna Marini.
        Manage time-related keywords.
        '''

        if not hasattr(self, 'location'):
            self.earthlocation()

        param = self.params['obstime']
        timekey = self.header[param]
        
        if 'MJD' in param:
            obstime = Time(timekey, format='mjd')
        elif param == 'JD':
            obstime = Time(timekey, format='jd')
        elif isinstance(param, list):
            obstime = Time(Time(timekey[0]).unix+timekey[1])
        elif 'DATE' in param:
            obstime = Time(timekey)
        else:
            obstime = Time(timekey)

        # For sidereal time
        obstime.location = self.location

        # For other methods
        self.obstime = obstime
        return obstime


    def skycoord(self):
        '''By Anna Marini.
        Manage keywords related to radec coordinates.
        '''

        if not hasattr(self, 'location'):
            self.earthlocation()

        if not hasattr(self, 'obstime'):
            self.time()

        unit = (u.hourangle, u.deg)

        if self.params['ra'] and self.params['dec']:
            ra  = self.header[self.params['ra']]
            dec = self.header[self.params['dec']]
            
            coord = SkyCoord(ra=ra, dec=dec, unit=unit)

        elif self.in_head('OBJECT'):
            obj = self.header['OBJECT']
            try:
                coord = SkyCoord.from_name(obj)
            except:
                print("Object not found in catalog")
                print("Provide ra dec or check object name")
                pass
        else:
            print("Not radec nor object found in header")
            pass

        # For alt, az
        coord.location = self.location
        coord.obstime = self.obstime

        # For other methods
        self.coord = coord
        return coord


    def meteo(self):
        '''
        Manage keywords related to weather conditions.
        TBD: Add to alt az params.
        '''

        self.temperature = None  # 20*u.Celsius
        self.humidity = None # 0-1
        self.pressure = None    # 1000*u.hpa
        self.wavelength = None  # 550*u.nm
        if self.in_head('XTEMP'):
            self.temperature = self.header['XTEMP']
        if self.in_head('HUMIDITY'):
            self.humidity = self.header['HUMIDITY']/100
        if self.in_head('ATMOSBAR'):
            self.pressure = self.header['ATMOSBAR']


    def detector(self):
        '''
        Manage keywords related to the detector.
        '''
                    
        # elif self.in_head('CCDSUM'):
        #     binning = list(map(int, self.header['CCDSUM'].split(' ')))

        if self.params["binning"]:
            binning = [self.header[self.params["binning"][0]],
                       self.header[self.params["binning"][1]]]
        else:
            binning = [1, 1]

        if self.params["scale"]:
            scale = self.params["scale"]
        else:
            scale = 1

        plate = scale * binning[0]
        self.plate = plate
        return plate
        
    def wcs(self):
        '''From Anna Marini.
        Provide WCS keywords to convert pixel coordinates of the
        files to sky coordinates. It uses the rotational matrix
        obtained in previous function (which_instrument)
        '''

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'plate'):
            self.detector()

        angle = 0
        flip = 1
        if self.header['INSTRUME'] == 'Mexman':
            angle = np.pi/2
            flip = -1

        cd = self.plate*np.array([[np.cos(angle), np.sin(angle)*flip],
                                  [np.sin(angle), np.cos(angle)]])

        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = cd
        w.wcs.crval = [self.coord.ra.deg,
                       self.coord.dec.deg]
        w.wcs.crpix = [self.header['NAXIS1']/2,
                       self.header['NAXIS2']/2]

        # For other methods
        self.w = w
        return w


    def in_head(self, s):
        '''
        Check if keywords in the header.
        '''
        if isinstance(s, list):
            h = [elem in self.header for elem in s ]
        else:
            h = s in self.header
        return h


    def sethead(self, key, val):
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
            comm = None
            form = '{}'
        else:
            form = '{'+card[0]["format"]+'}'
            comm = card[0]["comment"]
        
        if is_number(val):
            value = form.format(val)
            value = to_number(value)
        else:
            value = val
  
        print(key, val, "→", form, "→", value)
          
        if not card:
            self.header[key] = value
        else:
            self.header[key] = (value, card[0]["comment"])


    def newhead(self):

        if not hasattr(self, 'plate'):
            self.detector()

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'w'):
            self.wcs()

        # location
        self.sethead("latitude", self.location.lat.deg)
        self.sethead("longitud", self.location.lon.deg)
        self.sethead("altitude", int(self.location.height.to_value()) )

        # obstime
        self.sethead("mjd-obs", self.obstime.mjd)
        self.sethead("jd", self.obstime.jd)
        self.sethead("date-obs", self.obstime.isot[:-4])
        self.sethead("st", self.obstime.sidereal_time("mean").hourangle)
        self.sethead("equinox", self.obstime.jyear_str)
        
        # detector
        self.sethead("plate", self.plate)

        # coord
        #self.sethead("ra", self.coord.ra.to_string(unit="hourangle",sep=":"))
        #self.sethead("dec", self.coord.dec.to_string(sep=":"))
        self.sethead("ra", self.coord.ra.deg)
        self.sethead("dec", self.coord.dec.deg)

        self.sethead("alt", self.coord.altaz.alt.deg)
        self.sethead("az", self.coord.altaz.az.deg)
        self.sethead("airmass", self.coord.altaz.secz.value)
        self.sethead("zdist", self.coord.altaz.zen.deg)
        
        sun_radec = get_sun(self.obstime)
        sun_altaz = sun_radec.transform_to(self.coord.altaz)
        moon_radec = get_moon(self.obstime)
        moon_altaz = moon_radec.transform_to(self.coord.altaz)

        self.sethead("sunalt", sun_altaz.alt.deg)
        self.sethead("sundist", sun_radec.separation(self.coord).deg)
        self.sethead("moonalt", moon_altaz.alt.deg)
        self.sethead("moondist", moon_radec.separation(self.coord).deg)

        # wcs
        self.header.extend(self.w.to_header(), update=True)
        
        #self.header.extend(w.to_header(), update=True)
        # self.header.rename_keyword("PC1_1", "CD1_1", force=True)
