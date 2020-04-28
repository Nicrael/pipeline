#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import json

# Our modules
import reduction as r

class observatory():

    def __init__(self, **kwargs):
        '''
        Set default parameters.
        '''

        # Header keywords with defaults
        self.ra = kwargs.get('ra', 'RA')
        self.dec = kwargs.get('dec', 'DEC')
        self.obj = kwargs.get('obj', 'OBJECT')
        self.obstime = kwargs.get('obstime', 'MJD-OBS')

        # Header keywords optional
        self.binning = kwargs.get('binning')
        self.temperature = kwargs.get('temperature')
        self.pressure = kwargs.get('pressure')
        self.humidity = kwargs.get('humidity')

        # Constants and defaults
        self.lon =  kwargs.get('lon')
        self.lat =  kwargs.get('lat')
        self.alt =  kwargs.get('alt', 0)
        self.unit =  kwargs.get('unit', (u.hour, u.deg))
        self.scale =  kwargs.get('scale', 0.5) # in binning 1,

        # No header by default
        self.head = None

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
        self.head = r.get_fits_header(value)


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
            print("No (RA DEC) and no OBJECT in header")
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

        plate = self.scale * bins[0]

        self.bins = bins
        self.plate = plate
        return plate


    def wcs(self):
        '''By Anna Marini.
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


    def newhead(self):
        nh = fits.PrimaryHDU().header

        if not hasattr(self, 'plate'):
            self.detector()

        if not hasattr(self, 'coord'):
            self.skycoord()

        if not hasattr(self, 'w'):
            self.wcs()

        c = self.coord

        # location
        sethead(nh, "latitude", c.location.lat.deg)
        sethead(nh, "longitud", c.location.lon.deg)
        sethead(nh, "altitude", int(c.location.height.to_value()) )

        # obstime
        sethead(nh, "mjd-obs", c.obstime.mjd)
        sethead(nh, "jd", c.obstime.jd)
        sethead(nh, "date-obs", c.obstime.isot[:-4])
        sethead(nh, "lst", c.obstime.sidereal_time("mean").hour)
        sethead(nh, "equinox", c.obstime.jyear)

        # detector
        sethead(nh, "plate", self.plate)

        # coord
        #sethead(nh, "ra", c.ra.to_string(unit="hourangle",sep=":"))
        #sethead(nh, "dec", c.dec.to_string(sep=":"))
        sethead(nh, "ra", c.ra.deg)
        sethead(nh, "dec", c.dec.deg)

        sethead(nh, "alt", c.altaz.alt.deg)
        sethead(nh, "az", c.altaz.az.deg)
        sethead(nh, "airmass", c.altaz.secz.value)
        sethead(nh, "zdist", c.altaz.zen.deg)

        sun_radec = get_sun(c.obstime)
        sun_altaz = sun_radec.transform_to(c.altaz)
        moon_radec = get_moon(c.obstime)
        moon_altaz = moon_radec.transform_to(c.altaz)

        sethead(nh, "sunalt", sun_altaz.alt.deg)
        sethead(nh, "sundist", sun_radec.separation(c).deg)
        sethead(nh, "moonalt", moon_altaz.alt.deg)
        sethead(nh, "moondist", moon_radec.separation(c).deg)

        if hasattr(self, "_filename"):
            sethead(nh, "FULLPATH", self.filename) # .split("/")[-1])

        # wcs
        #nh.extend(self.w.to_header(), update=True)

        #self.header.extend(w.to_header(), update=True)
        # self.header.rename_keyword("PC1_1", "CD1_1", force=True)

        nh.add_comment("Created by "+type(self).__name__)

        self.nh = nh


def sethead(head, key, val):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file

    In [30]: 'Hello, {:.2f}'.format(123.129)
    Out[30]: 'Hello, 123.13'

    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-array.json') as jfile:
        head_format = json.load(jfile)# ['primary']

        #card = [fits.Card(**c) for c in head_format if key == c['keyword'] ][0]
        #form = '' if not card else card[0]["format"]

    for h in head_format :
        if key == h["keyword"]:
            form = h["value"]
            if 'd' in form:
                value = int(round(float(val)))
            elif 'f' in form:
                value = round(val, [ int(v) for v in h["value"] if v.isdigit() ][0])
            else: # 's' in form:
                value = f'{val:{h["value"]}}'
            print(key, val, "→", form, "→", value)
            head[key] = (value, h["comment"])
        # else:
        #     print(f'{key} not in dict')

    return head


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
        print(o.nh)


if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 3 :    # C'è anche lo [0] che è il nome del file :)
        print("Usage:  "+sys.argv[0]+" <instrument as in json> <list of FITS files>")
        sys.exit()

    main()
