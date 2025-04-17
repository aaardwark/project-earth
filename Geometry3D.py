import math
import numbers
from fractions import Fraction
"""
:# ETC - only validation code; bare bones
:# WIP - work in progress; involved math/logic
:# FIX - theory issue
:# BUG - code not working
"""

class Point:
    def __init__(self, x, y, z):
        self.x = Fraction(x)
        self.y = Fraction(y)
        self.z = Fraction(z)

    @staticmethod
    def from_vector(vector):
        if type(vector) != Vector:
            raise TypeError('Point.from_vector requires a Vector.')
        return Point(vector.x, vector.y, vector.z)

    def __repr__(self):
        return f'Point({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __str__(self):
        return f'({float(self.x)}, {float(self.y)}, {float(self.z)})'
    
    def __eq__(self, point):
        return type(point)==Point and self.x==point.x and self.y==point.y and self.z==point.z

    def __abs__(self):
        return float(math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z))
    
    def __add__(self, vector):
        if type(vector) != Vector:
            return NotImplemented
        return Point(self.x+vector.x, self.y+vector.y, self.z+vector.z)

    def __sub__(self, point):
        if type(point) != Point:
            return NotImplemented
        return Vector(self.x-point.x, self.y-point.y, self.z-point.z)
        
    def distance_to(self, geobj):
        if type(geobj) == Point:
            x_diff = self.x - geobj.x
            y_diff = self.y - geobj.y
            z_diff = self.z - geobj.z
            return math.sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff)
        elif type(geobj) == Line:
            NotImplemented
        elif type(geobj) == Plane:
            normal = geobj.perpendicular_line_through(self)
            point_on_plane = normal.intersect(geobj)
            return self.distance_to(point_on_plane)
        else:
            NotImplemented

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
    
    def __eq__(self, vector):
        return type(vector)==Vector and self.x==vector.x and self.y==vector.y and self.z==vector.z
    
    def __abs__(self):
        return float(math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z))
    
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
            raise TypeError('Argument to Vector.dot_product must be a vector')
        return self.x*vector.x + self.y*vector.y + self.z+vector.z
    
    def cross_product(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.cross_product must be a vector')
        return Vector(self.y*vector.z - self.z*vector.y, self.z*vector.x - self.x*vector.z, self.x*vector.y - self.y*vector.x)

    def angle_to(self, vector):
        if type(vector) != Vector:
            raise TypeError('Argument to Vector.angle_to must be a vector')
        elif not (abs(self) and abs(vector)):
            raise ValueError('Can only find the angle between two nonzero vectors.')
        cos_angle = self.dot_product(vector)/(abs(self) * abs(vector))
        return math.acos(cos_angle)


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
        
    @staticmethod
    def from_coefficients(x_coeff, y_coeff, z_coeff, const):
        if not (x_coeff and y_coeff and z_coeff):
            raise ValueError('Plane.from_coefficients requires at least one nonzero coefficient.')
        normal_vector = Vector(x_coeff, y_coeff, z_coeff)
        point = Point.from_vector( ORIGIN - const*normal_vector/(abs(normal_vector)**2) )
        return Plane(point, normal_vector)
    
    @staticmethod
    def from_points(point1, point2, point3):
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
    
    def __eq__(self, plane):
        if type(plane)==Plane:
            return (self.a/plane.a == self.b/plane.b == self.c/plane.c == self.d/plane.d)
        return False

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

        elif type(geobj) == Plane:# WIP
            cross_vector = Vector(self.a, self.b, self.c).cross_product(Vector(geobj.a, geobj.b, geobj.c))
            if abs(cross_vector) == 0:
                if self == geobj:
                    return self
                else:
                    return None
            else:
                
            
        
        elif type(geobj) == Sphere:
            NotImplemented
            
        else:
            NotImplemented

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
            raise TypeError('Line constructor requires one Point and one nonzero Vector argument.')
        self.point = point
        self.vector = vector
    
    @staticmethod
    def from_points(point1, point2):
        if type(point1)!=Point or type(point2)!=Point or point1==point2:
            raise TypeError('Line.from_points requires two different Point arguments.')
        return Line(point1, abs(point2 - point1))
    
    def __repr__(self):
        return f'Line({repr(self.point)}, {repr(self.vector)})'
        
    def __str__(self):
        sLp = self.point
        sLv = self.vector
        parametrise_x = f'{float(sLp.x)} + {float(sLv.x)}t'
        parametrise_y = f'{float(sLp.y)} + {float(sLv.y)}t'
        parametrise_z = f'{float(sLp.z)} + {float(sLv.z)}t'
        return f'({parametrise_x}, {parametrise_y}, {parametrise_z})'

    def __eq__(self, line):
        if type(line) == Line:
            return line.has_point(self.point) and self.vector.angle_to(line.vector)==0
        return False
        
    def has_point(self, point):
        if type(point) != Point:
            raise TypeError('Argument to Line.has_point must be a Point.')
        return self.vector.angle_to(point - self.point) == 0

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
            NotImplemented
        else:
            NotImplemented
    
    def project_onto_plane(self, plane):
        if type(plane) != Plane:
            raise TypeError('Argument to Line.project_onto_plane must be a Plane')
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
    
    def __eq__(self, sphere):
        if type(sphere) != Sphere:
            return False
        else:
            return (self.centre==sphere.centre and self.radius==sphere.radius)
    
    def has_point(self, point):
        if type(point) != Point:
            raise TypeError('Argument to Sphere.has_point must be a Point.')
        sSc = self.centre
        x_diff = point.x - sSc.x
        y_diff = point.y - sSc.y
        z_diff = point.z - sSc.z
        return ( (x_diff*x_diff + y_diff*y_diff + z_diff*z_diff) == (self.radius*self.radius) )

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
            NotImplemented
        elif type(geobj) == Sphere:
            NotImplemented
        else:
            NotImplemented

    def tangent_plane_at(self, point):
        if type(point)!=Point or not self.has_point(point):
            raise TypeError('Argument to Sphere.tangent_plane_at must be a Point lying on the sphere.')
        return Plane(point, self.centre - point)



def solve_quadratic(pwr2_coeff, pwr1_coeff, pwr0_coeff):
    if not (isinstance(pwr2_coeff, numbers.Real) and isinstance(pwr1_coeff, numbers.Real) and isinstance(pwr0_coeff, numbers.Real)):
        raise TypeError('All coefficients must be real numbers')
    
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