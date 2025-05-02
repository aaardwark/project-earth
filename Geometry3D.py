from math import *
import numbers
from fractions import Fraction
"""
:# ETC - only validation code; bare bones
:# WIP - work in progress; involved math/logic
:# FIX - theory issue
:# BUG - code not working
"""
"""
class Ray (and class Segment?)
class Spheroid? class Circle/Ellipse?
implement intersections
rotation transform of a point - any line, any centre
"""

class Point:
    def __init__(self, x, y, z):
        self.x = Fraction(x)
        self.y = Fraction(y)
        self.z = Fraction(z)

    @classmethod
    def from_vector(cls, vector):
        if type(vector) != Vector:
            raise TypeError('Point.from_vector requires a Vector.')
        return Point(vector.x, vector.y, vector.z)

    def __repr__(self):
        return f'Point({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __str__(self):
        return f'({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __eq__(self, other):
        return type(other)==Point and self.x==other.x and self.y==other.y and self.z==other.z

    def __abs__(self):
        return float(sqrt(self.x*self.x + self.y*self.y + self.z*self.z))
    
    """A Point can have another Point subtracted from it, giving a Vector.
    A Point can also have a Vector added to it, giving a Point."""
    def __add__(self, vector):
        if type(vector) != Vector:
            raise TypeError('Can only add Points to a Point.')
        return Point(self.x+vector.x, self.y+vector.y, self.z+vector.z)

    def __sub__(self, point):
        if type(point) != Point:
            raise TypeError('Can only subtract Points from a Point.')
        return Vector(self.x-point.x, self.y-point.y, self.z-point.z)
        
    def distance_to(self, geobj):
        if type(geobj) == Point:
            return abs(self - geobj)
        elif type(geobj) == Line:
            common_plane = Plane(self, geobj.vector)
            nearest_point = common_plane.intersect(geobj)
            return abs(self - nearest_point)
        elif type(geobj) == Plane:
            normal = geobj.perpendicular_line_through(self)
            point_on_plane = normal.intersect(geobj)
            return self.distance_to(point_on_plane)
        else:
            raise NotImplementedError(f'Cannot find distance from a Point to a {geobj.__class__.__name__}')

ORIGIN = Point(0,0,0)



class Vector:
    def __init__(self, x, y, z):
        self.x = Fraction(x)
        self.y = Fraction(y)
        self.z = Fraction(z)
    
    def __repr__(self):
        return f'Vector({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __str__(self):
        return f'({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __eq__(self, other):
        return type(other)==Vector and self.x==other.x and self.y==other.y and self.z==other.z
    
    def __abs__(self):
        return float(sqrt(self.x*self.x + self.y*self.y + self.z*self.z))
    
    def __add__(self, vector):
        if type(vector) != Vector:
            return NotImplemented
        return Vector(self.x+vector.x, self.y+vector.y, self.z+vector.z)

    def __sub__(self, vector):
        if type(vector) != Vector:
            return NotImplemented
        return Vector(self.x-vector.x, self.y-vector.y, self.z-vector.z)

    def __mul__(self, factor):
        if not isinstance(factor, numbers.Real):
            return NotImplemented
        return Vector(self.x*factor, self.y*factor, self.z*factor)
    
    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)
    
    def as_unit(self):
        return self * (1/abs(self))

    def dot_product(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.dot_product must be a Vector')
        return self.x*vector.x + self.y*vector.y + self.z+vector.z
    
    def cross_product(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.cross_product must be a Vector')
        return Vector(self.y*vector.z - self.z*vector.y, self.z*vector.x - self.x*vector.z, self.x*vector.y - self.y*vector.x)

    def angle_to(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.angle_to must be a Vector')
        elif not (abs(self) and abs(vector)):
            raise ValueError('Can only find the angle between two nonzero Vectors.')
        cos_angle = self.dot_product(vector)/(abs(self) * abs(vector))
        return acos(cos_angle)
    
    def is_parallel_to(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.is_parallel_to must be a Vector')
        return self.angle_to(vector) in {0,pi}


I = Vector(1,0,0)
J = Vector(0,1,0)
K = Vector(0,0,1)


class Plane:
    def __init__(self, point, normal_vector):
        if type(point)!=Point or type(normal_vector)!=Vector or abs(normal_vector)==0:
            raise TypeError('Plane constructor requires one Point and one nonzero Vector argument.')
        self.a = normal_vector.x
        self.b = normal_vector.y
        self.c = normal_vector.z
        self.d = 0 - normal_vector.x*point.x - normal_vector.y*point.y - normal_vector.z*point.z
        
    @classmethod
    def from_coefficients(cls, x_coeff, y_coeff, z_coeff, const):
        if not (x_coeff and y_coeff and z_coeff):
            raise ValueError('Plane.from_coefficients requires at least one nonzero coefficient.')
        normal_vector = Vector(x_coeff, y_coeff, z_coeff)
        point = Point.from_vector( ORIGIN - const*normal_vector/(abs(normal_vector)**2) )
        return Plane(point, normal_vector)
    
    @classmethod
    def from_points(cls, point1, point2, point3):
        if not type(point1) == type(point2) == type(point3) == Point:
            raise TypeError('Plane.from_points requires three distinct Points.')
        elif (point1 - point2).angle_to(point2 - point3) == 0:
            raise ValueError('Plane.from_points requires three non-collinear Points.')
        else:
            v1 = point3 - point1
            v2 = point3 - point2
            normal_vector = v1.cross_product(v2)
            return Plane(point3, normal_vector)

    def __repr__(self):
        return f'Plane.from_coefficients({self.a}, {self.b}, {self.c}, {self.d})'

    def __str__(self):
        return f'{float(self.a)}x + {float(self.b)}y + {float(self.c)}z + {float(self.d)} = 0'
    
    def __eq__(self, other):
        if type(other)==Plane:
            return self.normal_vector()*other.d == other.normal_vector()*self.d
        return False
    
    def normal_vector(self):
        return Vector(self.a, self.b, self.c)

    def is_parallel_to(self, plane):
        if type(plane) != Plane:
            raise TypeError(f'Argument to Plane.is_parallel_to must be a Plane')
        return self.normal_vector().is_parallel_to(plane.normal_vector())
        
    def has_point(self, point):
        if type(point) != Point:
            raise TypeError('Argument to Plane.has_point must be a Point.')
        return (self.a*point.x + self.b*point.y + self.c*point.z + self.d == 0)

    def intersect(self, geobj):
        if type(geobj) == Line:
            gLp = geobj.point
            gLv = geobj.vector
            parameter_coeff = self.a*gLv.x + self.b*gLv.y + self.c*gLv.z
            constant = self.a*gLp.x + self.b*gLp.y + self.c*gLp.z + self.d

            if parameter_coeff == 0:
                if constant == 0:
                    return geobj
                return None
            parameter_sol = -constant/parameter_coeff
            sol_x = gLv.x*parameter_sol + gLp.x
            sol_y = gLv.y*parameter_sol + gLp.y
            sol_z = gLv.z*parameter_sol + gLp.z
            return Point(sol_x, sol_y, sol_z)

        elif type(geobj) == Plane:
            cross_vector = self.normal_vector().cross_product(geobj.normal_vector())
            if abs(cross_vector) == 0:
                if self == geobj:
                    return self
                else:
                    return None
            else:
                raise NotImplementedError
            
        
        elif type(geobj) == Sphere:
            raise NotImplementedError
            
        else:
            raise NotImplementedError

    def perpendicular_line_through(self, point):
        if type(point)!=Point:
            raise TypeError('Argument to Plane.perpendicular_line_through must be a Point.')
        return Line(point, Vector(self.a, self.b, self.c))

YZ_PLANE = Plane(ORIGIN, I)
ZX_PLANE = Plane(ORIGIN, J)
XY_PLANE = Plane(ORIGIN, K)



class Line:
    def __init__(self, point, vector):
        if type(point)!=Point or type(vector)!=Vector or abs(vector)==0:
            raise TypeError(f'{self.__class__.__name__} constructor requires one Point and one nonzero Vector argument.')
        self.point = point
        self.vector = vector
    
    @classmethod
    def from_points(cls, point1, point2):
        if type(point1)!=Point or type(point2)!=Point or point1==point2:
            raise TypeError(f'{cls.__name__}.from_points requires two different Point arguments.')
        return cls(point1, abs(point2 - point1))
    
    def __repr__(self):
        return f'{self.__class__.__name__}({repr(self.point)}, {repr(self.vector)})'
        
    def __str__(self):
        sLp = self.point
        sLv = self.vector
        parametrise_x = f'{float(sLp.x)} + {float(sLv.x)}t'
        parametrise_y = f'{float(sLp.y)} + {float(sLv.y)}t'
        parametrise_z = f'{float(sLp.z)} + {float(sLv.z)}t'
        return f'({parametrise_x}, {parametrise_y}, {parametrise_z})'

    def __eq__(self, other):
        if type(other) == Line:
            return other.has_point(self.point) and self.vector.is_parallel_to(other.vector)
        return False

    def is_parallel_to(self, line):
        if type(line) != Line:
            raise TypeError(f'Argument to {self.__class__.__name__}.is_parallel_to must be a Line')
        return self.vector.is_parallel_to(line.vector)
        
    def has_point(self, point):
        if type(point) != Point:
            raise TypeError(f'Argument to {self.__class__.__name__}.has_point must be a Point.')
        return self.vector.is_parallel_to(point - self.point)

    def intersect(self, geobj):
        sLp = self.point
        sLv = self.vector
        if type(geobj) == Line:
            if self == geobj:
                return self
            gLp = geobj.point
            gLv = geobj.vector
            does_intersect = (sLv.z*gLv.x - sLv.x*gLv.z)(sLv.x*gLp.y - sLv.x*sLp.y - sLv.z*gLp.x + sLv.z*sLp.x) == (sLv.y*gLv.x - sLv.x*gLv.y)(sLv.x*gLp.z - sLv.x*sLp.z - sLv.z*gLp.x + sLv.z*sLp.x)

            if does_intersect:
                s = (sLv.x*gLp.y - sLv.x*sLp.y - sLv.z*gLp.x + sLv.z*sLp.x)/(sLv.y*gLv.x - sLv.x*gLv.y)
                return Point(gLv.x*s + gLp.x, gLv.y*s + gLp.y, gLv.z*s + gLp.z)
            return None

        elif type(geobj) == Plane:
            return geobj.intersect(self)
        elif type(geobj) == Sphere:
            raise NotImplementedError
        else:
            raise NotImplementedError
    
    def project_onto_plane(self, plane):
        if type(plane) != Plane:
            raise TypeError(f'Argument to {self.__class__.__name__}.project_onto_plane must be a Plane')
        sLp = self.point
        sLv = self.vector
        point_coeff = (plane.a*sLp.x + plane.b*sLp.y + plane.c*sLp.z + plane.d)/(plane.a**2 + plane.b**2 + plane.c**2)
        point = Point(sLp.x - plane.a*point_coeff, sLp.y - plane.b*point_coeff, sLp.z - plane.c*point_coeff)
        vector_coeff = (plane.a*sLv.x + plane.b*sLv.y + plane.c*sLv.z)/(plane.a**2 + plane.b**2 + plane.c**2)
        vector = Vector(sLv.x - plane.a*vector_coeff, sLv.y - plane.b*vector_coeff, sLv.z - plane.c*vector_coeff)
        if abs(vector) == 0:
            return point
        else:
            return Line(point, vector)

X_AXIS = Line(ORIGIN, I)
Y_AXIS = Line(ORIGIN, J)
Z_AXIS = Line(ORIGIN, K)


class Ray(Line):
    # inherit __init__
    # inherit classmethod from_points
    # inherit __repr__

    def __str__(self):
        return super().__str__(self) + ' {t ≥ 0}'
    
    def __eq__(self, other):
        if type(other) == Ray:
            return other.point==self.point and self.vector.angle_to(other.vector)==0
        return False
        

class Sphere:
    # generalise to Spheroid
    def __init__(self, centre, radius):
        if type(centre) != Point:
            raise TypeError('First argument to Sphere constructor must be a Point.')
        self.centre = centre
        self.radius = abs(Fraction(radius))

    def __repr__(self):
        return f'Sphere({repr(self.centre)}, {float(self.radius)})'

    def __str__(self):
        sSc = self.centre
        return f'(x - {float(sSc.x)})^2 + (y - {float(sSc.y)})^2 + (z - {float(sSc.z)})^2 = {float(self.radius*self.radius)}'
    
    def __eq__(self, other):
        if type(other) == Sphere:
            return (self.centre==other.centre and self.radius==other.radius)
        return False

    def volume(self):
        return Fraction(pi * self.radius**3 * 4/3)
    
    def surface_area(self):
        return Fraction(pi * self.radius**2 * 4)
    
    def has_point(self, point):
        if type(point) != Point:
            raise TypeError('Argument to Sphere.has_point must be a Point.')
        sSc = self.centre
        x_diff = point.x - sSc.x
        y_diff = point.y - sSc.y
        z_diff = point.z - sSc.z
        return (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff) == (self.radius*self.radius)

    def intersect(self, geobj):
        sSc = self.centre
        if type(geobj) == Line:
            gLp = geobj.point
            gLv = geobj.vector
            t_sqr_coeff = gLv.x**2 + gLv.y**2 + gLv.z**2
            const = (gLp.x - sSc.x)**2 + (gLp.y - sSc.y)**2 + (gLp.z - sSc.z)**2 - self.radius**2

            t_coeff = 2 * gLv.x * (gLp.x - sSc.x) 
            t_coeff += 2 * gLv.y * (gLp.y - sSc.y) 
            t_coeff += 2 * gLv.z * (gLp.z - sSc.z) 

            return solve_quadratic(t_sqr_coeff, t_coeff, const)

        elif type(geobj) == Plane:
            raise NotImplementedError
        elif type(geobj) == Sphere:
            raise NotImplementedError
        else:
            raise NotImplementedError

    def tangent_plane_at(self, point):
        if type(point)!=Point or not self.has_point(point):
            raise TypeError('Argument to Sphere.tangent_plane_at must be a Point lying on the sphere.')
        return Plane(point, self.centre - point)



def solve_quadratic(pwr2_coeff, pwr1_coeff, pwr0_coeff):
    if not (isinstance(pwr2_coeff, numbers.Real) and isinstance(pwr1_coeff, numbers.Real) and isinstance(pwr0_coeff, numbers.Real)):
        raise TypeError('All coefficients in a quadratic must be real numbers')
    
    a = Fraction(pwr2_coeff)
    b = Fraction(pwr1_coeff)
    c = Fraction(pwr0_coeff)
    
    discriminant = (b*b) - (4*a*c)
    if discriminant<0:
        return None
    elif discriminant>0:
        return ( (-b + discriminant)/(2*a) , (-b - discriminant)/(2*a) )
    else:
        return ( (-b)/(2*a), )

rad = lambda angle_in_degrees: angle_in_degrees*pi/180
deg = lambda angle_in_radians: angle_in_radians*180/pi

def arctrig(cosv, sinv):
        asinset = {round(asin(sinv), 8), round((pi - asin(sinv))%pi, 8)}
        acosset = {round(acos(cosv), 8), round(-acos(cosv), 8)}
        return Fraction( asinset.intersection(acosset).pop() )


def rotation_xy(point, centre, angle):
    if type(point) != Point:
        raise TypeError('Can only rotate Points.')
    elif type(centre) != Point:
        raise TypeError('Centre of rotation has to be a Point.')
    
    translated_pt = point - (centre + Vector(0,0,point.z))
    modulus = abs(translated_pt)
    start_angle = arctrig(translated_pt.x/modulus, translated_pt.y/modulus)
    rotated_angle = start_angle + angle
    rotated_pt = Point( modulus*cos(rotated_angle), modulus*sin(rotated_angle) )
    return rotated_pt + centre

def rotate_about_axis(point, axis_name, angle):
    if type(point) != Point:
        raise TypeError('Can only rotate Points.')

    if axis_name == 'x':
        pt_projected = (point.y,point.z)
    elif axis_name == 'y':
        pt_projected = (point.z,point.x)
    elif axis_name == 'z':
        pt_projected = (point.x,point.y)
    else:
        raise NameError(f'No axis called: {axis_name}')
    
    modulus = sqrt(pt_projected[0]**2 + pt_projected[1]**2)
    start_angle = arctrig(pt_projected.x/modulus, pt_projected.y/modulus)
    rotated_angle = start_angle + angle
    pt_rotated = ( modulus*cos(rotated_angle), modulus*sin(rotated_angle) )

    if axis_name == 'x':
        return Point(point.x, pt_rotated[0], pt_rotated[1])
    elif axis_name == 'y':
        return Point(pt_rotated[1], point.y, pt_rotated[0])
    else:
        return Point(pt_rotated[0], pt_rotated[1], point.z)
