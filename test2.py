# System modules
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
import astropy.units as u


# h = fits.getheader("examples/201610080221o.fits.fz", 1)
# l = EarthLocation(lon= -115.4869,
#                   lat= 31.0292,
#                   height= 2790
# )

# t = Time(h['jd'],
#          format="jd",
#          location=l)

# coor = SkyCoord(ra=h['ra'],
#                 dec=h['dec'],
#                 unit=(u.hourangle, u.deg),
#                 location=l,
#                 obstime=t)


h = fits.getheader("examples/DFOSC_FASU.2010-05-10T04:41:33.000.fits", 0)
l = EarthLocation(lon="-70:44:14.7926",
                  lat= "-29:15:27.8194",
                  height=2375
)

t = Time(h['mjd-obs'],
         format="mjd",
         location=l)

coor = SkyCoord(ra=h['ra'],
                dec=h['dec'],
                unit=(u.hourangle, u.deg),
                location=l,
                obstime=t)
