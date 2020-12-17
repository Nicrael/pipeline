#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Complement the FITS header with useful information.
'''

# System modules
from astropy import log
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.coordinates import FK5
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import json
import astropy.units as u
import numpy as np

# Local modules
from fits import get_fits_header
from naming import hist


class Observatory():
    '''
    Observatory class.
    '''

    def __init__(self, **kwargs):
        log.info("Set default parameters: {kwargs}")

        # Header keywords with defaults
        self.ra = kwargs.get('ra', 'RA')
        self.dec = kwargs.get('dec', 'DEC')
        self.obj = kwargs.get('obj', 'OBJECT')
        self.obstime = kwargs.get('obstime', 'JD')

        # Header keywords optional
        self.binning = kwargs.get('binning')
        self.temperature = kwargs.get('temperature')
        self.pressure = kwargs.get('pressure')
        self.humidity = kwargs.get('humidity')

        # Constants and defaults
        self.lon = kwargs.get('lon')
        self.lat = kwargs.get('lat')
        self.alt = kwargs.get('alt', 0)
        self.unit = kwargs.get('unit', (u.hour, u.deg))
        self.scale = kwargs.get('scale', 0.5)  # in binning 1,

        # No header by default
        self.head = fits.PrimaryHDU().header


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
        '''
        Property for filename
        '''
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        '''
        Setter for filename
        '''
        self._filename = value
        self.head = get_fits_header(value)


    def skycoord(self):
        '''
        Manage and fix sky coordinates
        '''
        head = self.head.copy()

        # time
        if self.obstime not in head:
            head['JD'] = Time.now().jd

        # earthlocation
        location = EarthLocation(lon=self.lon,
                                 lat=self.lat,
                                 height=self.alt)

        timekey = head[self.obstime]
        # log.warning(timekey)
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

        obstime.location = location

        # skycoord
        if self.ra in head and self.dec in head:
            coord = SkyCoord(ra=head[self.ra],
                             dec=head[self.dec],
                             location=location,
                             obstime=obstime,
                             equinox=obstime.jyear_str[:-2],
                             # frame='icrs',
                             frame='fk5',
                             unit=self.unit)

        elif self.obj in head:
            coord = SkyCoord.from_name(head[self.obj],
                                       #location=location,
                                       #obstime=obstime,
                                       # frame='icrs',
                                       #equinox=obstime.jyear_str[:-2]
                                       frame='fk5')

        else:
            log.error("No (RA DEC) and no OBJECT in header")
            return

        coord = coord.transform_to(FK5(equinox='J2000'))

        # Meteo
        if self.temperature and self.temperature in head:
            coord.temperature = head[self.temperature] * u.deg_C
        if self.pressure and self.pressure in head:
            coord.pressure = head[self.pressure] * u.hPa
        if self.humidity and self.humidity in head:
            coord.relative_humidity = head[self.humidity] / 100

        # coord.obstime=obstime

        self.coord = coord
        return coord

    def detector(self):
        '''
        Manage keywords related to the detector.
        '''
        head = self.head.copy()

        if self.binning and self.binning[0] in head:
            bins = [head[self.binning[0]],
                    head[self.binning[1]]]
        else:
            bins = [1, 1]

        plate = (self.scale * bins[0])*u.arcsec.to(u.deg)

        self.bins = bins
        self.plate = plate
        return plate

    def wcs(self):
        '''By Anna Marini.
        Provide WCS keywords to convert pixel coordinates of the
        files to sky coordinates. It uses the rotational matrix
        obtained in previous function (which_instrument)
        '''
        head = self.head.copy()

        # if not hasattr(self, 'coord'):
        self.coord = self.skycoord()

        # if not hasattr(self, 'plate'):
        self.plate = self.detector()

        angle = 0
        flip = 1
        if 'INSTRUME' in head and head['INSTRUME'] == 'Mexman':
            log.warning("Found Mexman! Rotating WCS!")
            angle = np.pi/2
            flip = 1

        wcss = WCS(head)
        wcss.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        wcss.wcs.cd = self.plate*np.array([[np.cos(angle), flip*np.sin(angle)],
                                           [np.sin(angle), np.cos(angle)]])

        wcss.wcs.crval = [self.coord.ra.deg,
                          self.coord.dec.deg]
        wcss.wcs.crpix = [(head['NAXIS1']+1)/2,
                          (head['NAXIS2']+1)/2]

        # For other methods
        self.wcss = wcss
        return wcss

    def newhead(self, header=False):
        '''
        Build the new header with all useful information.
        '''
        
        if header:
            self.head = header

        nhd = self.head.copy()

        # if not hasattr(self, 'plate'):
        self.detector()

        # if not hasattr(self, 'coord'):
        self.skycoord()

        # if not hasattr(self, 'wcss'):
        self.wcs()

        coord = self.coord

        # location
        nhd["longitud"] = coord.location.lon.deg
        nhd["latitude"] = coord.location.lat.deg
        nhd["altitude"] = int(coord.location.height.to_value())

        # # location (hierarch test)
        # nhd["TEL GEOLON"] = coord.location.lon.deg
        # nhd["TEL GEOLAT"] = coord.location.lat.deg
        # nhd["TEL GEOELEV"] = int(coord.location.height.to_value())

        # obstime
        nhd["mjd-obs"] = coord.obstime.mjd
        #nh["jd"] = coord.obstime.jd
        nhd["date-obs"] = coord.obstime.isot  # [:-4]

        midnight = coord.obstime.iso.split()[0]
        nhd["utc"] = coord.obstime.unix - Time(midnight).unix

        # sid = coord.obstime.sidereal_time("mean").hour # MUST NOT BE J2000
        #nhd["lst"] = sid*u.hour.to(u.s)

        nhd["equinox"] = coord.equinox.jyear  # should be 2000 for fk5

        # detector
        nhd["plate"] = self.plate

        if self.obj in nhd:
            nhd[self.obj] = nhd[self.obj]

        # coord
        #nhd["ra"] = coord.ra.to_string(unit="hourangle",sep=":")
        #nhd["dec"] = coord.decoord.to_string(sep=":")
        nhd["radesys"] = coord.frame.name.upper()
        nhd["ra"] = coord.ra.deg
        nhd["dec"] = coord.dec.deg
        # nh["ha"] = sid - coord.ra # MUST NOT BE IN J2000

        nhd["alt"] = coord.altaz.alt.deg
        nhd["az"] = coord.altaz.az.deg
        nhd["airmass"] = coord.altaz.secz.value
        nhd["zdist"] = coord.altaz.zen.deg

        sun_radec = get_sun(coord.obstime)
        sun_altaz = sun_radec.transform_to(coord.altaz)
        moon_radec = get_moon(coord.obstime)
        moon_altaz = moon_radec.transform_to(coord.altaz)

        nhd["sunalt"] = sun_altaz.alt.deg
        nhd["sundist"] = sun_radec.separation(coord).deg
        nhd["moonalt"] = moon_altaz.alt.deg
        nhd["moondist"] = moon_radec.separation(coord).deg

        #nhd.extend(w.to_header(), update=True)

        # if not hasattr(self, 'w'):
        #     log.warn("calling wcs!!!")
        #     #if "NAXIS1" in nh:
        #     self.wcs()
        #     nhd.extend(self.w.to_header(), update=True)

        if hasattr(self, "_filename"):
            nhd["FULLPATH"] = self.filename  # .split("/")[-1])

        nhd.add_history(hist(__name__))

        nhd = sethead(nhd)

        self.nhd = nhd
        return nhd


def solver(pattern, ra=False, dec=False, scale=False):
    '''
    Calls the solve-field command from the astrometry.net
    debian package.
    '''

    import os

    cmd = f'solve-field {pattern} \
    --crpix-center \
    --downsample 2 \
    --no-plots \
    --dir solved \
    --overwrite'

    if (ra and dec):
        cmd += f' \
        --radius 0.1 \
        --ra {ra} \
        --dec {dec}'

    if scale:
        cmd += f' \
        --scale-units arcsecperpix \
        --scale-low {float(scale)*0.9} \
        --scale-high {float(scale)*1.1}'

    print(cmd)

    os.system(cmd)


def sethead(head):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file
    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-test.json') as jfile:
        header_dict = json.load(jfile)  # ['primary']

    template = fits.Header([tuple(d) for d in header_dict if d[0] in head])

    bastard_keywords = {"COMMENT", "HISTORY"}

    for i in template:
        #log.info(f"formatting {n}")
        if i not in bastard_keywords:
            form = template[i]
            val = head[i]
            value = formatter(val, form)
            template.set(*(i, value, ))

    for i in bastard_keywords:
        if i in template:
            value = head[i]
            for j in value:
                template.add_history(j)

    diff = fits.HeaderDiff(head, template)
    if diff:
        log.warning(f"Not in dictionary:{diff.diff_keywords}")

    return template


def formatter(val, form):
    '''
    Take a value and format it according to
    the format string taken by the header dictionary.
    '''
    #log.info(f"formatting {val} to {form}")

    if 'd' in form:  # integer
        value = int(round(float(val)))
    elif 'f' in form:  # float
        decimals = [int(v) for v in form if v.isdigit()]
        decimals = [0] if not decimals else decimals
        value = round(val, decimals[0])
    elif 'b' in form:  # boolean
        value = bool(val)
    elif 's' in form:  # string
        value = f'{val:{form}}'
    else:  # 'x' in form: ## history or comment
        value = val
        #log.info(f"Other: {value}")
    return value


def init_observatory(instrument="Mexman"):
    '''
    Init observatory using a json config file.
    '''

    log.info(f'Loading {instrument}')

    with open('./instruments.json') as jfile:
        instruments = json.load(jfile)

    instrument = instruments[instrument]
    return instrument  # s, observatory(**instrument)
