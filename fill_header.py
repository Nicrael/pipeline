#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %load_ext autoreload
# %autoreload 2

# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import get_sun, get_moon
from astropy.time import Time
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
        
    
    def config(self, filename=None):
        if filename is not None:
            self.filename = filename
        self.coordinates()
        self.location()
        self.times()
        self.detector()
        self.meteo()
        self.altaz()

        
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


    # def target(self):
    #     if self.in_head('OBJECT'):
    #         target = self.header['OBJECT']
    #     else:
    #         target = None
    #     return SkyCoord.from_name(target)


    def coordinates(self):
        value=self.params
        if not value:
            pass # List is empty
        else:
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

        return skycoord


    def location(self):
        value = self.params
        if not value:
            pass # List is empty
        else:
            if all(self.in_head(['LATITUDE','LONGITUD','ALTITUDE'])):
                lat = self.header[value['location'][0]]
                lon = self.header[value['location'][1]]
                alt = self.header[value['location'][2]]
            else:
                lat = value['location'][0]
                lon = value['location'][1]
                alt = value['location'][2]

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
        value = self.params
        if not value:
            pass # List is empty
        else:
            timekey = self.header[value['obstime']]
            if 'MJD' in value['obstime']:
                obstime = Time(timekey, format='mjd')
            elif value['obstime'] == 'JD':
                obstime = Time(timekey, format='jd')
            elif 'DATE' in value['obstime']:
                obstime = Time(timekey)
            elif isinstance(value['obstime'], list):
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
        
        target_radec = self.coordinates()
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

            
    def in_head(self, s):
        '''
        Check if a keyword is in the header.
        '''
        if isinstance(s, list):
            h = [elem in self.header for elem in s ]
        else:
            h = s in self.header
        return h    


    
        
    # # json array
    # [item for item in j if item.get('id')=='Mexman'

    # # json object
    # j['Mexman']




# Mettere tutto il JSON in un oggetto Python per poter trattare i dati come credi


# def fill_header(filename):

#     observing_location = EarthLocation(lat=41.3*u.deg,
#                                        lon=-74*u.deg,
#                                        height=100*u.m)

#     observing_time = Time('2017-02-05 20:12:18')
#     aa = AltAz(location=observing_location, obstime=observing_time)

#     obj_radec = SkyCoord('4h42m', '-38d6m50.8s')
#     obj_altaz = object_radec.transform_to(aa)

#     temperature = None # 20*u.Celsius
#     pressure = None    # 1000*u.hpa
#     humidity = None    # 0.2
#     wavelength = None  # 550*u.nm

#     zdist = object_altaz.zen # only in dfosc
#     airmass = object_altaz.secz # only in mexman

#     sun_radec = get_sun(observing_time)
#     sun_altaz = sun_radec.transform_to(aa)

#     moon_radec = get_moon(observing_time)
#     moon_altaz = moon_radec.transform_to(aa)


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


