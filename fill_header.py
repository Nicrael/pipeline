#!/usr/bin/env python3
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
        '''By Davide Ricci.
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

        # if all(self.in_head(['LONGITUD','LATITUDE','ALTITUDE'])):
        #     lon = self.header[param[0]]
        #     lat = self.header[param[1]]
        #     alt = self.header[param[2]]
        # else:
        
        lon = param[0]
        lat = param[1]
        alt = param[2]

        # if self.header['INSTRUME'] == 'Mexman':
        #     if isinstance(self.header['LONGITUD'], str):
        #         lon = "-"+lon
                
        location = EarthLocation(lat=lat, lon=lon, height=alt)

        #if (self.in_head('OBSERVAT'):
        #    earthlocation = EarthLocation.of_site(self.header('OBSERVAT'))
        #   earthlocation = EarthLocation.of_address("")

        # For header
        self.sethead("latitude", location.lat.deg)
        self.sethead("longitud", location.lon.deg)
        self.sethead("altitude", int(location.height.to_value()) )

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

        # if self.in_head('EQUINOX'):
        #     equitime = self.header['EQUINOX']
        #     if str(equitime).startswith('J'):
        #         equinox = Time(equitime)
        #     else:
        #         equinox = Time(equitime, format='jyear')
        # else:
        #     equinox = obstime

        # For sidereal time
        obstime.location = self.location

        # For header
        self.sethead("mjd-obs", obstime.mjd)
        self.sethead("jd", obstime.jd)
        self.sethead("date-obs", obstime.isot[:-4])
        self.sethead("st", obstime.sidereal_time("mean").hourangle)
        self.sethead("equinox", obstime.jyear_str)

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


        if self.params['ra'] and self.params['dec']:
            ra  = self.header[self.params['ra']]
            dec = self.header[self.params['dec']]

            if self.header['INSTRUME'] == 'DFOSC_FASU':
                #example dfosc: 14.32572 decimal hours
                coord = SkyCoord(ra=ra, dec=dec,
                                 unit=(u.hourangle, u.deg))
                
            elif is_number(ra):
                coord = SkyCoord(ra=ra, dec=dec,
                                 unit=(u.deg, u.deg))
            else:
                #example mexman: 18:56:10.8
                coord = SkyCoord(ra=ra, dec=dec,
                                 unit=(u.hourangle, u.deg))

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

        # if self.in_head('ZDIST'):
        #     zdist = self.header['ZDIST'] # example: dfosc
        # else:
        #    zdist = coord.altaz.zen.deg

        # if self.in_head('AIRMASS'):
        #     airmass = self.header['AIRMASS'] # example: mexman
        # else:
        #     airmass = coord.altaz.secz.value

        sun_radec = get_sun(self.obstime)
        sun_altaz = sun_radec.transform_to(coord.altaz)
        moon_radec = get_moon(self.obstime)
        moon_altaz = moon_radec.transform_to(coord.altaz)

        # For header
        #self.sethead("ra", coord.ra.to_string(unit="hourangle",sep=":"))
        #self.sethead("dec", coord.dec.to_string(sep=":"))
        self.sethead("ra", coord.ra.deg)
        self.sethead("dec", coord.dec.deg)

        self.sethead("alt", coord.altaz.alt.deg)
        self.sethead("az", coord.altaz.az.deg)
        self.sethead("airmass", coord.altaz.secz.value)
        self.sethead("zdist", coord.altaz.zen.deg)

        self.sethead("sunalt", sun_altaz.alt.deg)
        self.sethead("sundist", sun_radec.separation(coord).deg)
        self.sethead("moonalt", moon_altaz.alt.deg)
        self.sethead("moondist", moon_radec.separation(coord).deg)

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

        if all(self.in_head(['CCDXBIN','CCDYBIN'])):
            binning = [self.header['CCDXBIN'],
                       self.header['CCDXBIN']]
        elif self.in_head('CCDSUM'):
            binning = list(map(int, self.header['CCDSUM'].split(' ')))

        elif all(self.in_head(['XBINNING','YBINNING'])):
            binning = [self.header['XBINNING'],
                       self.header['YBINNING']]
        else:
            binning = [1, 1]

        # For header
        self.sethead("xbin", binning[0])
        self.sethead("ybin", binning[1])
        if self.params["scale"] is not None:
            self.sethead("scale", self.params["scale"])
        else:
            self.sethead("scale", 1)


    def wcs(self):
        '''From Anna Marini.
        Provide WCS keywords to convert pixel coordinates of the
        files to sky coordinates. It uses the rotational matrix
        obtained in previous function (which_instrument)

        '''

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'scale'):
            self.detector()

        plate = (self.header["scale"] * self.header["xbin"])/3600
        angle = 0
        flip = 1
        if self.header['INSTRUME'] == 'Mexman':
            angle = np.pi/2
            flip = -1

        cd = np.array([[plate*np.cos(angle), plate*np.sin(angle)*flip],
                       [plate*np.sin(angle), plate*np.cos(angle)]])

        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = cd
        w.wcs.crval = [self.coord.ra.deg,
                       self.coord.dec.deg]
        w.wcs.crpix = [self.header['NAXIS1']/2,
                       self.header['NAXIS2']/2]
        #o, in alternativa, x_target e y_target date in input

        # For header
        self.header.extend(w.to_header(), update=True)
        # self.header.rename_keyword("PC1_1", "CD1_1", force=True)
        # self.header.rename_keyword("PC1_2", "CD1_2", force=True)
        # self.header.rename_keyword("PC2_1", "CD2_1", force=True)
        # self.header.rename_keyword("PC2_2", "CD2_2", force=True)

        # For other methods
        self.w = w
        return w


    def testarray(self):
        '''
        Just a test.
        '''
        with open('cerbero-merged-array.json') as cm:
            ccc = json.load(cm)

        aaa = fits.PrimaryHDU()
        #if any([var.startswith('%') for var in ccc['primary'][key]['comment'].split()]): # variable in comment
        for item in ccc['primary']:
            val = item['default'] if 'default' in item else None
            sss = fits.Card(item["name"], val, item["comment"])
            aaa.header.extend([sss], update=True)

        sss = [ ("OARPAF "+c["name"], c["default"], c["comment"]) for c in ccc['hierarch'] if 'default' in c ]
        aaa.header.extend(sss, update=True)


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
            
