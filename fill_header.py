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
from reduction import get_fits_header, is_number


class observatory():

    def __init__(self, filename=None):
        '''
        Set default parameters.
        '''
        # No default filename, header or instrument
        self._filename = None
        self.original = None # original header
        self.header = None
        self.instrument = None
        self.params = None

        with open('instruments.json') as json_file:
            self.instruments = json.load(json_file)
            # # json array
            # [item for item in j if item.get('id')=='Mexman'
            # # json object
            # j['Mexman']

        if filename is not None:
            self.filename = filename

    def config(self, filename=None):
        '''By Anna Marini.
        Runs all methods to have all parameters.
        '''
        if filename is not None:
            self.filename = filename
        self.location()
        self.times()
        self.coordinates()
        self.altaz()
        self.detector()
        self.meteo()
        self.wcs()
        
        #self.header.extend(w.to_header(), update=True)

    @property
    def filename(self):
        '''Laboratory image file name'''
        return self._filename

    @filename.setter # On new file, update data
    def filename(self, value):
        self.original =  get_fits_header(value)
        self.header = self.original.copy()
        self._filename = value
        if self.in_head('INSTRUME'):
            self.instrument = self.header['INSTRUME']
        else:
            self.instrument = 'default'
        self.params = self.instruments[self.instrument]
        self.exptime = self.header[self.params['exptime']]


    def coordinates(self):
        '''By Anna Marini.
        Manage keywords related to radec coordinates.
        '''
        
        if self.params['ra'] and self.params['dec']:
            ra  = self.header[self.params['ra']]
            dec = self.header[self.params['dec']]

            if is_number(ra):
                #example dfosc: 14.32572
                skycoord = SkyCoord(ra=ra, dec=dec,
                                    unit=(u.deg, u.deg))
            else:
                #example mexman: 18:56:10.8
                skycoord = SkyCoord(ra=ra, dec=dec,
                                    unit=(u.hourangle, u.deg))
        elif self.in_head('OBJECT'):
            target = self.header['OBJECT']
            try:
                skycoord = SkyCoord.from_name(target)
            except:
                print("Object not found in catalog")
                print("Provide ra dec or check object name")
                exit()
        else:
            print("Not radec nor object found in header")
            exit()

        # For header
        #self.header["ra"] = skycoord.ra.to_string(unit="hourangle",sep=":")
        #self.header["dec"] = skycoord.dec.to_string(sep=":")
        self.header["ra"] = skycoord.ra.deg
        self.header["dec"] = skycoord.dec.deg

        # For other methods
        self.skycoord = skycoord
        return skycoord


    def location(self):
        '''By Davide Ricci.
        Manage keywords related to local parameters.
        '''
        param = self.params['location']
        if all(self.in_head(['LONGITUD','LATITUDE','ALTITUDE'])):
            lon = self.header[param[0]]
            lat = self.header[param[1]]
            alt = self.header[param[2]]
        else:
            lon = param[0]
            lat = param[1]
            alt = param[2]
            
        if self.instrument == 'Mexman':
            if isinstance(self.header['LONGITUD'], str):
                lon = "-"+lon
            
        earthlocation = EarthLocation(lat=lat, lon=lon, height=alt)

        #if (self.in_head('OBSERVAT'):
        #    earthlocation = EarthLocation.of_site(self.header('OBSERVAT'))
        #   earthlocation = EarthLocation.of_address("")

        # For header
        self.header["latitude"] = earthlocation.lat.deg
        self.header["longitud"] = earthlocation.lon.deg
        self.header["altitude"] = int(earthlocation.height.to_value())

        # For other methods
        self.earthlocation = earthlocation
        return earthlocation


    def times(self):
        '''By Anna Marini.
        Manage time-related keywords.
        '''

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

        if self.in_head('EQUINOX'): #check format! Must contain "J"
            equitime = self.header['EQUINOX']
            if str(equitime).startswith('J'):
                equinox = Time(equitime)
            else:
                equinox = Time(equitime, format='jyear')
        else:
            equinox = obstime

        # For header
        self.header["mjd-obs"] = obstime.mjd
        self.header["jd"] = obstime.jd
        self.header["date-obs"] = obstime.isot[:-4]
        self.header["equinox"] = equinox.jyear_str

        # For other methods
        self.obstime = obstime
        return obstime


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
        self.header["xbin"] = binning[0]
        self.header["ybin"] = binning[1]
        if self.params["scale"] is not None:
            self.header["scale"] = self.params["scale"]
        else:
            self.header["scale"]  = 1


    def meteo(self):
        '''
        Manage keywords related to weather conditions.
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


    def altaz(self):
        '''
        Manage keywords related to local geographic position.
        '''

        if not hasattr(self, 'skycoord'):
            self.coordinates()
        if not hasattr(self, 'earthlocation'):
            self.location()
        if not hasattr(self, 'obstime'):
            self.times()

        target_radec = self.skycoord
        observing_location = self.earthlocation
        observing_time = self.obstime

        generic_altaz = AltAz(location=observing_location,
                              obstime=observing_time)
        #add meteo stuff for altaz

        target_altaz = target_radec.transform_to(generic_altaz)

        if self.in_head('ZDIST'):
            zdist = self.header['ZDIST'] # example: dfosc
        else:
            zdist = target_altaz.zen.deg

        if self.in_head('AIRMASS'):
            airmass = self.header['AIRMASS'] # example: mexman
        else:
            airmass = target_altaz.secz.value

        sun_radec = get_sun(observing_time)
        sun_altaz = sun_radec.transform_to(generic_altaz)
        moon_radec = get_moon(observing_time)
        moon_altaz = moon_radec.transform_to(generic_altaz)

        # For header
        self.header["alt"] = target_altaz.alt.deg
        self.header["az"] = target_altaz.az.deg
        self.header["airmass"] = airmass
        self.header["zdist"] = zdist
        self.header["sunalt"] = sun_altaz.alt.deg
        self.header["sundist"] = sun_radec.separation(target_radec).deg
        self.header["moonalt"] = moon_altaz.alt.deg
        self.header["moondist"] = moon_radec.separation(target_radec).deg

        # For other methods
        self.generic_altaz = generic_altaz
        return generic_altaz


    def wcs(self):
        '''From Anna Marini.
        Provide WCS keywords to convert pixel coordinates of the
        files to sky coordinates. It uses the rotational matrix
        obtained in previous function (which_instrument)

        '''

        if not hasattr(self, 'skycoord'):
            self.coordinates()

        if not hasattr(self, 'scale'):
            self.detector()

        plate = (self.header["scale"] * self.header["xbin"])/3600
        angle = 0
        flip = 1
        if self.instrument == 'Mexman':
            angle = np.pi/2
            flip = -1

        cd = np.array([[plate*np.cos(angle), plate*np.sin(angle)*flip],
                       [plate*np.sin(angle), plate*np.cos(angle)]])

        w = WCS()
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = cd
        w.wcs.crval = [self.skycoord.ra.deg,
                       self.skycoord.dec.deg]
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


    def in_head(self, s):
        '''
        Check if a keyword is in the header.
        '''
        if isinstance(s, list):
            h = [elem in self.header for elem in s ]
        else:
            h = s in self.header
        return h


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


def main():
    '''
    Main function.
    '''
    pattern = sys.argv[1:] # "1:" stands for "From 1 on".
    for filename in pattern:
        print(filename)
        o=observatory()
        o.config(filename)

if __name__ == '__main__':
    '''
    If called as a script.
    '''
    import sys

    if len(sys.argv) < 2 :
        print("Usage:  "+sys.argv[0]+" <list of fits files>")
        sys.exit()

    main()
