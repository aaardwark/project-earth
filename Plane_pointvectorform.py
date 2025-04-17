
class Plane:
    def __init__(self, point, normal_vector):
        if type(point)!=Point or type(normal_vector)!=Vector or abs(normal_vector)==0:
            raise TypeError('Plane constructor requires one Point and one nonzero Vector argument.')
        self.point = point
        self.normal_vector = normal_vector
        self.d = 0 - normal_vector.x*point.x - normal_vector.y*point.y - normal_vector.z*point.z
        
    @staticmethod
    def from_coefficients(x_coeff, y_coeff, z_coeff, const):
        if not (x_coeff and y_coeff and z_coeff):
            raise ValueError('Plane constructor requires at least one nonzero coefficient.')
        normal_vector = Vector(x_coeff, y_coeff, z_coeff)
        point = Point(ORIGIN - normal_vector.unit_vector()*const)
        return Plane(point, normal_vector)
    
    @staticmethod
    def from_points(point1, point2, point3):
        if points_are_collinear(point1, point2, point3):
            raise ValueError('Cannot construct a Plane out of collinear points.')
        else:
            p1, p2, p3 = point1, point2, point3
            x_coeff = (p1.z-p3.z)*(p2.y-p3.y) - (p2.z-p3.z)*(p1.y-p3.y)
            y_coeff = (p2.z-p3.z)*(p1.x-p3.x) - (p1.z-p3.z)*(p2.x-p3.x)
            z_coeff = (p1.y-p3.y)*(p2.x-p3.x) - (p1.x-p3.x)*(p2.y-p3.y)
            constant = 0 - (x_coeff)*p1.x - (y_coeff)*p1.y - (z_coeff)*p1.z

            return Plane.from_coefficients(x_coeff, y_coeff, z_coeff, constant)

    def __repr__(self):
        return f'Plane({repr(self.point)}, {repr(self.normal_vector)})'

    def __str__(self):
        sPv = self.normal_vector
        return f'{float(sPv.x)}x + {float(sPv.y)}y + {float(sPv.z)}z + {float(self.d)} = 0'
    
    def __eq__(self, plane):
        if type(plane)==Plane:
            return (self.normal_vector.unit_vector() == plane.normal_vector.unit_vector() and self.has_point(plane.point))
        return False

    def has_point(self, point):
        if type(point) != Point:
            raise TypeError('Argument to Plane.has_point must be a Point.')
        sPp = self.point
        sPv = self.normal_vector
        return (sPv.x*point.x + sPv.y*point.y + sPv.z*point.z == sPv.x*sPp.x + sPv.y*sPp.y + sPv.z*sPp.z)

    def intersect(self, geobj):
        sPv = self.normal_vector
        if type(geobj) == Line:
            gLp = geobj.point
            gLv = geobj.vector
            parameter_coeff = sPv.x*gLv.x + sPv.y*gLv.y + sPv.z*gLv.z
            constant = sPv.x*gLp.x + sPv.y*gLp.y + sPv.z*gLp.z + self.d

            if parameter_coeff == 0:
                if constant == 0:
                    return geobj
                return None
            parameter_sol = -constant/parameter_coeff
            sol_x = geobj.vector.x*parameter_sol + geobj.point.x
            sol_y = geobj.vector.y*parameter_sol + geobj.point.y
            sol_z = geobj.vector.z*parameter_sol + geobj.point.z
            return Point(sol_x, sol_y, sol_z)

        elif type(geobj) == Plane:# WIP
            NotImplemented
        
        elif type(geobj) == Sphere:
            NotImplemented
            
        else:
            NotImplemented

    def perpendicular_line_through(self, point):
        if type(point)!=Point:
            raise TypeError('Argument to Plane.perpendicular_line_through must be a Point.')
        return Line(point, self.normal_vector)
