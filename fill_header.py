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
        head = self.head

        # earthlocation
        location = EarthLocation(lon=self.lon,
                                 lat=self.lat,
                                 height=self.alt)

        # time
        if self.obstime not in self.head:
            self.head['JD'] = Time.now().jd

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
                             frame='fk5',
                             unit=self.unit)
        elif self.obj in head:
            coord = SkyCoord.from_name(head[self.obj], frame='fk5')
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
            bins = [head[self.binning[0]],
                    head[self.binning[1]]]
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

        w = WCS(head)
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

        nh = self.head

        if not hasattr(self, 'plate'):
            self.detector()

        if not hasattr(self, 'coord'):
            self.skycoord()

        c = self.coord

        # location
        nh["longitud"] = c.location.lon.deg
        nh["latitude"] = c.location.lat.deg
        nh["altitude"] = int(c.location.height.to_value())

        # location (hierarch test)
        nh["tel geolon"] = c.location.lon.deg
        nh["tel geolat"] = c.location.lat.deg
        nh["tel geoelev"] = int(c.location.height.to_value())

        # obstime
        nh["mjd-obs"] = c.obstime.mjd
        nh["jd"] = c.obstime.jd
        nh["date-obs"] = c.obstime.isot[:-4]

        midnight = c.obstime.iso.split()[0]

        nh["utc"] = c.obstime.unix - Time(midnight).unix
        nh["lst"] = c.obstime.sidereal_time("mean").hour

        nh["equinox"] = c.equinox.jyear # should be 2000 for fk5

        # detector
        nh["plate"] = self.plate

        if self.obj in nh:
            nh[self.obj] = nh[self.obj]

        # coord
        #nh["ra"] = c.ra.to_string(unit="hourangle",sep=":")
        #nh["dec"] = c.dec.to_string(sep=":")
        nh["radecsys"] = c.frame.name
        nh["ra"] = c.ra.deg
        nh["dec"] = c.dec.deg

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

        if not hasattr(self, 'w'):
            if "NAXIS1" in nh:
                self.wcs()
                nh.extend(self.w.to_header(), update=True)

        if hasattr(self, "_filename"):
            nh["FULLPATH"] = self.filename # .split("/")[-1])

        nh.add_history("Created by "+self.newhead.__qualname__)

        nh = sethead2(nh)

        self.nh = nh


def sethead2(head):
    ''' By Davide Ricci.
    Add a keyword in the header with comment and format
    taken from the json config file
    In [47]: f'Hello, {123.129:.2f}'
    Out[51]: 'Hello, 123.13'
    '''

    with open('cerbero-merged-test.json') as jfile:
        header_dict = json.load(jfile)# ['primary']

    nh = fits.Header([tuple(d) for d in header_dict if d[0] in head ] )

    for n in nh:
        form = nh[n]
        val = head[n]
        #comm = nh.get_comment(n)
        if 'd' in form:
            value = int(round(float(val)))
        elif 'f' in form:
            decimals = [ int(v) for v in form if v.isdigit() ]
            decimals = [0] if not decimals else decimals
            value = round(val, decimals[0])
        elif 'b' in form:
            value = bool(val)
        else: # 's' in form:
            value = f'{val:{form}}'
        #print(n, val, "→", form, "→", value)
        nh[n] = (value, )

    return nh


# def sethead(head):
#     ''' By Davide Ricci.
#     Add a keyword in the header with comment and format
#     taken from the json config file
#     In [47]: f'Hello, {123.129:.2f}'
#     Out[51]: 'Hello, 123.13'
#     '''

#     with open('cerbero-merged-dict.json') as jfile:
#         header_dict = json.load(jfile)# ['primary']

#     #card = [fits.Card(**c) for c in head_dict if key == c['keyword'] ][0]
#     #form = '' if not card else card[0]["format"]

#     for key in head :
#         key =  key.lower()
#         if key in header_dict: # keyword in the dictionary
#             val = head[key]
#             form = header_dict[key][0]  # Format in the dictionary
#             comm =  header_dict[key][1] # Comment in the dictionary
#             if 'd' in form:
#                 value = int(round(float(val)))
#             elif 'f' in form:
#                 decimals = [ int(v) for v in form if v.isdigit() ]
#                 value = round(val, decimals[0])
#             else: # 's' in form:
#                 value = f'{val:{form}}'
#             print(key, val, "→", form, "→", value)
#             head[key] = (value, comm)
#         else:
#             print(f'{key} not in dict')

#     return head



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
