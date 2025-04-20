from math import *
from fractions import Fraction
from Geometry3D import *
# ensure Geometry3D is accessible as a module

##### FRACTIONISE #####
piF = Fraction(pi)

sinF = lambda x: Fraction(sin(x))
cosF = lambda x: Fraction(cos(x))

sqrtF = lambda x: Fraction(sqrt(x))

##### ASSIGN CONSTANTS #####
earth = Sphere(ORIGIN, 1)
earth_radius = 6378 
# earth equatorial radius, in km to nearest km
# greater than polar radius

axial_tilt = Fraction(0.4091) 
# earth axial tilt, in radians to 4sf
earth_axis = Line(ORIGIN, Vector(1,0,tan(pi/2 - axial_tilt)))

eccentricity = Fraction(0.0167)
# eccentricity of earth-sun orbit, to 3sf
# ratio of focus-centre distance to semimajor length
period = Fraction(365.256)
# earth-sun orbital period, in days to nearest 1/1000 of a day
apoapsis = 15_210_000_000
periapsis = 14_710_000_000
# max and min earth-sun distace, in km

long_periapsis = Fraction(1.80)
# longitude of periapsis, ~103°, in radians to 3sf
# angle as seen at sun, between earth at sep equinox and at periapsis


date_input = None
time_input = None


##### FUNCTION DEFINITIONS #####
def arctrig(sinv, cosv):
    asinset = {round(asin(sinv), 8), round((piF - asin(sinv))%piF, 8)}
    acosset = {round(acos(cosv), 8), round(-acos(cosv), 8)}
    return Fraction( asinset.intersection(acosset).pop() )

def compute_eccentric_anomaly(mean_anomaly):
    # mean anomaly is equivalent to season_angle
    Ei = [0, 1] # some arbitrary values -piF< <piF
    while abs(Ei[1] - Ei[0]) > 0.0001:
        Ei.append(mean_anomaly + eccentricity*sin(Ei.pop(0)))
    return Ei[1]

def sun_position(days_from_periapsis):
    ### everything here is 2D, in the xy plane
    mean_anomaly = 2 * piF * days_from_periapsis / period
    eccentric_anomaly = compute_eccentric_anomaly(mean_anomaly)
    # angle relative to periapsis, coordinates relative to ellipse centre

    cartpos_rel_peraps = semimajor*(cosF(eccentric_anomaly) - eccentricity), semiminor*sinF(eccentric_anomaly)
    polpos_rel_peraps = arctrig(*cartpos_rel_peraps), sqrtF(cartpos_rel_peraps[0]**2 + cartpos_rel_peraps[1]**2)
    # cartesian and polar coordinates relative to central body, periapsis as baseline

    polpos_rel_nsolst = long_periapsis + piF/2 + polpos_rel_peraps[0], polpos_rel_peraps[1]
    cartpos_rel_nsolst = polpos_rel_nsolst[1]*cosF(polpos_rel_nsolst[0]), polpos_rel_nsolst[1]*sinF(polpos_rel_nsolst[1])
    # polar and cartesian coordinates relative to central body, northern solstice as baseline
    return cartpos_rel_nsolst

def days_since_periapsis(date):
    pass

sun = Point(*sun_position(days_since_periapsis(date_input)) , 0)

def corrected_longitude(time):
    pass

def point_at_latlong(latitude_deg, season_angle, time_angle):
    lat = latitude_deg * piF/180
    alpha = Fraction(axial_tilt)
    phi = Fraction(season_angle)
    T = Fraction(time_angle)

    latitude_plane = earth_axis.perpendicular_plane_through( Point(0, 0, sin(lat)/cos(alpha)) )
    latcircle_projection_xrad = Line( Point(cos(lat)*cos(-alpha), 0, cos(lat)*sin(-alpha)), earth_axis.vector).intersect(XY_PLANE).x

    long = compute_longitude_angle(phi, T)

    point_latcircle_projection = Point( cos(long)*latcircle_projection_xrad , sin(long) , 0 )
    point_latlong = Line(point_latcircle_projection, earth_axis.vector).intersect(latitude_plane)
    line_to_sun = Line.from_points(point_latlong, sun_position(phi))
    tangent_plane = earth.tangent_plane_at(point_latlong)



