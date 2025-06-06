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
NORTH_POLE = Point(0,0,EARTH_RADIUS)


def sun_position(mean_anomaly):
    "mean anomaly: angle representing the duration of the orbit, starting from periapsis, that has elapsed. returns a Point"
    # compute eccentric anomaly
    Ei = [0, 1] # arbitrary starting values; -pi< <pi
    while abs(Ei[1] - Ei[0]) > 0.0001:
        Ei.append( mean_anomaly + ECCENTRICITY*sin(Ei.pop(0)) )
    eccentric_anomaly = Ei.pop()

    position_1 = Point(SEMIMAJOR*(cos(eccentric_anomaly) - ECCENTRICITY), SEMIMINOR*sin(eccentric_anomaly), 0)
    "in ellipse-axes frame of reference, treating xy plane as orbital plane"

    position_2 = rotate_about_axis(position_1, 'z', radians(LONGITUDE_OF_PERIAPSIS-90))
    "in solstice/equinox frame of reference, treating xy plane as orbital plane"

    position_3 = rotate_about_axis(position_2, 'y', EARTH_AXIAL_TILT)
    "in solstice/equinox frame of reference, on actual orbital plane"

    return position_3


def location_on_earth(time, latitude_radians, longitude_ref):
    """time of day as the angle of rotation completed by the earth since 0000 hrs; radians \n
    latitude; degrees \n
    the sun - a reference point to define where on earth noon/midnight are; Point"""
    ref_longitude = -arctrig(longitude_ref.x, longitude_ref.y)
    corrected_longitude = ref_longitude - time
    longitude_vector = Vector(cos(corrected_longitude), sin(corrected_longitude), 0)
    
    axis_point_at_lat = Point(0,0, EARTH_RADIUS * sin(latitude_radians) )
    ray_through_location = Ray(axis_point_at_lat, longitude_vector)
    return ray_through_location.intersect(EARTH)


def compute_observed_angles(location_on_earth, position_of_sun):
    location_tangent = EARTH.tangent_plane_at(location_on_earth)

    north_ray = Ray.from_points(location_on_earth, NORTH_POLE)
    north_ray_proj = north_ray.project_onto_plane(location_tangent)
    sun_ray = Ray.from_points(location_on_earth, position_of_sun)
    sun_ray_proj = sun_ray.project_onto_plane(location_tangent)
    # typerror ray from equal points if location is given as n pole
    # have to handle case where sun is directly above location - projection will be a point

    bearing = north_ray_proj.vector.angle_to(sun_ray_proj.vector)
    # angle_to cannot differentiate symmetry


def parse_timestr_24hrs(timestr_24hrs):
    def invalid():
        raise
    truncated_time_str = timestr_24hrs[:4]
    invalid() if len(truncated_time_str) < 4 else 0
    try:
        time_int = int(invalid)
    except:
        invalid()
    invalid() if time_int < 0 or time_int > 2400 else 0

    hours, minutes = divmod(time_int, 60)
    return {'hours':hours, 'minutes':minutes}


def get_inputs():
    "should return time angle, latitude and mean anomaly as a tuple. rewrite according to how input is to be taken"
    time_inp = input('Enter the time in 24 hour format (only numbers - no spaces, dots, etc): ')
    latitude_inp = input('Enter your latitude in degrees: ')
    days_since_periapsis_inp = input('Enter number of days since periapsis: ')

    # process inputs and convert to appropriate format
    try:
        time_angle = None # TBD
        latitude_rad = radians(float(latitude_inp))
        mean_anomaly = (2*pi/PERIOD) * (int(days_since_periapsis_inp) % PERIOD)
        return time_angle, latitude_rad, mean_anomaly
    except Exception as e:
        print(f'{e.__class__.__name__}: {e.args[0]}')
        return get_inputs()
        

time_angular, latitude_radians, mean_anomaly = get_inputs()
SUN = sun_position()
LOCATION = location_on_earth()