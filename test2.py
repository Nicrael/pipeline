# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.coordinates import Longitude, Latitude,  Angle
from astropy.io import fits
from astropy.time import Time
import astropy.units as u

head = fits.getheader("examples/DFOSC_FASU.2010-05-10T04:41:33.000.fits", 0)

ra = 'ra'
dec = 'dec'
unit = (u.hour, u.deg)
obstime ='mjd-obs'

l = EarthLocation(lon="-70:44:14.7926",
                  lat= "-29:15:27.8194",
                  height=2375
)

t = Time(h['mjd-obs'],
         format="mjd",
         location=l)

coor = SkyCoord(ra=ra,
                dec=dec,
                unit=unit,
                location=l,
                obstime=t)
