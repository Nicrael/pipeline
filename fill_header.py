#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy import log
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.coordinates import FK5
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
import astropy.units as u
import numpy as np
import json

# Local modules
from fits import get_fits_header
from naming import hist


class observatory():

    def __init__(self, **kwargs):
        log.info(f"Set default parameters: {kwargs}")

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
        self.lon =  kwargs.get('lon')
        self.lat =  kwargs.get('lat')
        self.alt =  kwargs.get('alt', 0)
        self.unit =  kwargs.get('unit', (u.hour, u.deg))
        self.scale =  kwargs.get('scale', 0.5) # in binning 1,

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
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        self._filename = value
        self.head = r.get_fits_header(value)

    def skycoord(self):
        head = self.head.copy()

        # time
        if self.obstime not in head:
            head['JD'] = Time.now().jd

        # earthlocation
        location = EarthLocation(lon=self.lon,
                                 lat=self.lat,
                                 height=self.alt)

        timekey = head[self.obstime]
        #log.warning(timekey)
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
                             #frame='icrs',
                             frame='fk5',
                             equinox=obstime.jyear_str[:-2],
                             unit=self.unit)

        elif self.obj in head:
            coord = SkyCoord.from_name(head[self.obj],
                                       location=location,
                                       obstime=obstime,
                                       #frame='icrs',
                                       frame='fk5',
                                       equinox=obstime.jyear_str[:-2])

        else:
            log.error("No (RA DEC) and no OBJECT in header")
            return

        coord = coord.transform_to(FK5(equinox='J2000'))


        # Meteo
        if self.temperature and self.temperature in head:
            coord.temperature = head[self.temperature] *u.deg_C
        if self.pressure and self.pressure in head:
            coord.pressure = head[self.pressure] *u.hPa
        if self.humidity and self.humidity in head:
            coord.relative_humidity = head[self.humidity]/ 100

        #coord.obstime=obstime

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

        #if not hasattr(self, 'coord'):
        self.coord = self.skycoord()

        #if not hasattr(self, 'plate'):
        self.plate = self.detector()

        angle = 0
        flip = 1
        if 'INSTRUME' in head and head['INSTRUME'] == 'Mexman':
            log.warning("Found Mexman! Rotating WCS!")
            angle = np.pi/2
            flip = 1

        cd = self.plate*np.array([[np.cos(angle), flip*np.sin(angle)],
                                  [np.sin(angle), np.cos(angle)]])

        w = WCS(head)
        w.wcs.ctype = ["RA---TAN-SIP", "DEC--TAN-SIP"]
        w.wcs.cd = cd
        w.wcs.crval = [self.coord.ra.deg,
                       self.coord.dec.deg]
        w.wcs.crpix = [(head['NAXIS1']+1)/2,
                       (head['NAXIS2']+1)/2]

        # For other methods
        self.w = w
        return w


    def newhead(self, header=False):

        if header:
            self.head = header

        nh = self.head.copy()

        #if not hasattr(self, 'plate'):
        self.detector()

        #if not hasattr(self, 'coord'):
        self.skycoord()
        
        #if not hasattr(self, 'w'):
        w = self.wcs()

        c = self.coord

        # location
        nh["longitud"] = c.location.lon.deg
        nh["latitude"] = c.location.lat.deg
        nh["altitude"] = int(c.location.height.to_value())

        # # location (hierarch test)
        # nh["TEL GEOLON"] = c.location.lon.deg
        # nh["TEL GEOLAT"] = c.location.lat.deg
        # nh["TEL GEOELEV"] = int(c.location.height.to_value())

        # obstime
        nh["mjd-obs"] = c.obstime.mjd
        #nh["jd"] = c.obstime.jd
        nh["date-obs"] = c.obstime.isot #[:-4]

        midnight = c.obstime.iso.split()[0]
        nh["utc"] = c.obstime.unix - Time(midnight).unix

        #sid = c.obstime.sidereal_time("mean").hour # MUST NOT BE J2000
        #nh["lst"] = sid*u.hour.to(u.s)

        nh["equinox"] = c.equinox.jyear # should be 2000 for fk5

        # detector
        nh["plate"] = self.plate

        if self.obj in nh:
            nh[self.obj] = nh[self.obj]

        # coord
        #nh["ra"] = c.ra.to_string(unit="hourangle",sep=":")
        #nh["dec"] = c.dec.to_string(sep=":")
        nh["radesys"] = c.frame.name.upper()
        nh["ra"] = c.ra.deg
        nh["dec"] = c.dec.deg
        #nh["ha"] = sid - c.ra # MUST NOT BE IN J2000

        nh["alt"] = c.altaz.alt.deg
        nh["az"] = c.altaz.az.deg
        nh["airmass"] = c.altaz.secz.value
        nh["zdist"] = c.altaz.zen.deg

        sun_radec = get_sun(c.obstime)
        sun_altaz = sun_radec.transform_to(c.altaz)
        moon_radec = get_moon(c.obstime)
        moon_altaz = moon_radec.transform_to(c.altaz)

        nh["sunalt"] = sun_altaz.alt.deg
        nh["sundist"] = sun_radec.separation(c).deg
        nh["moonalt"] = moon_altaz.alt.deg
        nh["moondist"] = moon_radec.separation(c).deg

        #nh.extend(w.to_header(), update=True)

        # if not hasattr(self, 'w'):
        #     log.warn("calling wcs!!!")
        #     #if "NAXIS1" in nh:
        #     self.wcs()
        #     nh.extend(self.w.to_header(), update=True)

        if hasattr(self, "_filename"):
            nh["FULLPATH"] = self.filename # .split("/")[-1])

        nh.add_history(hist(__name__))

        nh = sethead(nh)

        self.nh = nh
        return nh
        

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
        header_dict = json.load(jfile)# ['primary']

    template = fits.Header([tuple(d) for d in header_dict if d[0] in head ] )

    bastard_keywords = {"COMMENT", "HISTORY"}

    for n in template:
        #log.info(f"formatting {n}")
        if n not in bastard_keywords:
            form = template[n]
            val = head[n]
            value = formatter(val, form)
            template.set(*(n, value, ))

    for bk in bastard_keywords:
        if bk in template:
            value = head[bk]
            for v in value:
                template.add_history(v)

    diff = fits.HeaderDiff(head, template)
    if diff:
        log.warn(f"Not in dictionary:{diff.diff_keywords}")

    return template


def formatter(val, form):
    '''
    Take a value and format it according to
    the format string taken by the header dictionary.
    '''
    #log.info(f"formatting {val} to {form}")


    if 'd' in form: # integer
        value = int(round(float(val)))
    elif 'f' in form: # float
        decimals = [ int(v) for v in form if v.isdigit() ]
        decimals = [0] if not decimals else decimals
        value = round(val, decimals[0])
    elif 'b' in form: # boolean
        value = bool(val)
    elif 's' in form: # string
        value = f'{val:{form}}'
    else:# 'x' in form: ## history or comment
        value = val
        #log.info(f"Other: {value}")
    return value


def init_observatory(instrument):
    log.info(f'Loading {instrument}')

    with open('./instruments.json') as jfile:
        instruments = json.load(jfile)

    instrument = instruments[instrument]
    return instruments, observatory(**instrument)
