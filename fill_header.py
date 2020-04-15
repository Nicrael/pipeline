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

    def __init__(self,filename=None):
        '''
        Set default parameters.
        '''
        # No default filename, header or instrument
        self._filename = None
        self.header = None
        self.instrument = None
        self.params = None
        
        with open('instruments.json') as json_file:
            self.instruments = json.load(json_file)
            # # json array
            # [item for item in j if item.get('id')=='Mexman'            
            # # json object
            # j['Mexman']
        
    
    def config(self, filename=None):
        if filename is not None:
            self.filename = filename
        self.coordinates()
        self.location()
        self.times()
        self.detector()
        self.meteo()
        self.altaz()
        self.wcs()
        #self.header.extend(w.to_header(), update=True)

        
    @property
    def filename(self):
        '''Laboratory image file name'''
        return self._filename
    
    @filename.setter # On new file, update data
    def filename(self, value):
        self.header = get_fits_header(value)
        self._filename = value
        if self.in_head('INSTRUME'):
            self.instrument = self.header['INSTRUME']
        else:
            self.instrument = 'default'
        self.params = self.instruments[self.instrument]
        self.exptime = self.header[self.params['exptime']]


    def coordinates(self):
        value = self.params
        if value['ra'] and value['dec']:
            ra  = self.header[value['ra']]
            dec = self.header[value['dec']]

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
            skycoord = SkyCoord.from_name(target)
        else:
            #example oarpaf: None
            skycoord = None
                
        # for header
        self.ra = skycoord.ra.to_string(unit="hourangle",sep=":")
        self.dec = skycoord.dec.to_string(sep=":")
        self.radeg = skycoord.ra.deg
        self.decdeg = skycoord.dec.deg

        self.skycoord = skycoord 
        return skycoord


    def location(self):
        param = self.params['location']
        if all(self.in_head(['LATITUDE','LONGITUD','ALTITUDE'])):
            lat = self.header[param[0]]
            lon = self.header[param[1]]
            alt = self.header[param[2]]
        else:
            lat = param[0]
            lon = param[1]
            alt = param[2]

        earthlocation = EarthLocation(lat=lat,
                                      lon=lon,
                                      height=alt)

        self.lat = earthlocation.lat.to_string(sep=":")
        self.lon = earthlocation.lon.to_string(sep=":")
        self.altitude = int(earthlocation.height.to_value())
            
        #if (self.in_head('OBSERVAT'):
        #    earthlocation = EarthLocation.of_site(self.header('OBSERVAT'))
        #   earthlocation = EarthLocation.of_address("")
           
        return earthlocation


    def times(self):
        param = self.params['obstime']
        timekey = self.header[param]
        if 'MJD' in value:
            obstime = Time(timekey, format='mjd')
        elif param == 'JD':
            obstime = Time(timekey, format='jd')
        elif 'DATE' in param:
            obstime = Time(timekey)
        elif isinstance(param, list):
            obstime = Time(Time(timekey[0]).unix+timekey[1])
        else:
            pass
            
        # for header
        self.mjd = obstime.mjd 
        self.jd = obstime.jd
        self.dateobs = obstime.isot
           
        return obstime


    def detector(self):

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

        # for header
        self.binning = binning
        self.scale = self.params["scale"] if self.params["scale"] is not None else 1 
        
    def meteo(self):
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
        observing_location = self.location()
        observing_time = self.times()
        generic_altaz = AltAz(location=observing_location,
                              obstime=observing_time)
        #add meteo stuff for altaz
        
        target_radec = self.skycoord
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

        # for header
        self.alt = target_altaz.alt.deg
        self.az = target_altaz.az.deg
        self.airmass = airmass
        self.zdist = zdist
        self.sunalt = sun_altaz.alt.deg
        self.sundist = sun_radec.separation(target_radec).deg
        self.moonalt = moon_altaz.alt.deg 
        self.moondist = moon_radec.separation(target_radec).deg
        
        return generic_altaz

    def wcs(self):
        '''From Anna Marini.
        Provides WCS keywords to convert pixel coordinates of the files
        to sky coordinates. It uses the rotational matrix obtained in previous
        function (which_instrument)'''

        plate = (self.scale * self.binning[0])/3600
        angle = 0
        flip = 1
        if self.instrument == 'Mexman':
            angle = np.pi/2
            flip = -1
            
        cd = np.array([[plate*np.cos(angle), plate*np.sin(angle)*flip],
                       [plate*np.sin(angle), plate*np.cos(angle)]])


        w = WCS(self.header)

        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.cd = cd

        w.wcs.crval = [self.skycoord.ra.deg,
                       self.skycoord.dec.deg]
        w.wcs.crpix = [self.header['NAXIS1']/2,
                       self.header['NAXIS2']/2]
        #o, in alternativa, x_target e y_target date in input

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

    

def main():
    '''
    Main function
    '''
    pattern = sys.argv[1:] # "1:" stands for "From 1 on".
    for filename in pattern:
        o=observatory()
        o.config(filename)
        print(o.__dict__)

if __name__ == '__main__':
    '''
    If called as a script
    '''
    import sys

    if len(sys.argv) < 2 :
        print("Usage:  "+sys.argv[0]+" <Parameters>")
        sys.exit()

    main()


