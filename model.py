from math import *
from fractions import Fraction
from Geometry3D import *
# ensure Geometry3D is accessible as a module

"https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html"
"J2000"
##### CONSTANTS #####

EARTH_RADIUS = 6378 
"earth equatorial radius; kilometres"

EARTH_AXIAL_TILT = Fraction(23.44) 
"earth axial tilt relative to orbital plane; 4sf; degrees"

APOAPSIS = 152_100_000
"maximum earth-sun distance; 10^5; km;"

PERIAPSIS = 147_100_000
"minimum earth-sun distance; 10^5; km"

SEMIMAJOR = 149_598_000
"semimajor axis of earth-sun orbit; 10^3; km"

SEMIMINOR = 149_577_000
"semiminor axis of earth-sun orbit; 10^3; km"

LONGITUDE_OF_PERIAPSIS = 103
"longitude of periapsis - angle as seen at sun, between earth at sep equinox and at periapsis; 3sf; degrees"

ECCENTRICITY = 0.0167
"eccentricity of earth-sun orbit - ratio of focus-centre distance to semimajor length; 3sf"

PERIOD = 365.256
"earth-sun orbital period; nearest 1/1000; days"

##### SETUP #####
EARTH = Sphere(ORIGIN, EARTH_RADIUS)
EARTH_AXIS = Z_AXIS
POINT_ON_EARTH = None

ORBITAL_PLANE = Plane.from_coefficients( tan(rad(EARTH_AXIAL_TILT)) , 0, -1, 0)

def compute_eccentric_anomaly(mean_anomaly):
    """mean anomaly is angle representing the duration of the orbit, starting from periapsis, that has been traversed \n
    eccentric anomaly is the angle at the centre of the orbit between periapsis and the projected point on the auxiliary circle"""
    Ei = [0, 1] # arbitrary starting values; -pi< <pi
    while abs(Ei[1] - Ei[0]) > 0.0001:
        Ei.append( mean_anomaly + ECCENTRICITY*sin(Ei.pop(0)) )
    return Ei.pop()

def sun_position(days_since_periapsis):
    "returns a Point"
    eccentric_anomaly = compute_eccentric_anomaly(2*pi*days_since_periapsis/PERIOD)

    position_1 = Point(SEMIMAJOR*(cos(eccentric_anomaly) - ECCENTRICITY), SEMIMINOR*sin(eccentric_anomaly), 0)
    "in ellipse-axes frame of reference, treating xy plane as orbital plane"

    position_2 = rotate_about_axis(position_1, 'z', rad(LONGITUDE_OF_PERIAPSIS-90))
    "in solstice/equinox frame of reference, treating xy plane as orbital plane"

    position_3 = rotate_about_axis(position_2, 'y', EARTH_AXIAL_TILT)
    "in solstice/equinox frame of reference, on actual orbital plane"

    return position_3


def point_on_earth(time, latitude, longitude_ref):
    """time of day as the angle of rotation completed by the earth since 0000 hrs; radians \n
    latitude; degrees \n
    the sun - a reference point to define where on earth noon/midnight are; Point"""
    latitude_plane = Plane.from_coefficients(0,0,-1, EARTH_RADIUS * sin(rad(latitude)))
    axis_point_at_lat = Point(0,0, EARTH_RADIUS * sin(rad(latitude)) )





days_since_periapsis = int(input('Enter number of days since periapsis: '))
SUN = sun_position(days_since_periapsis)
