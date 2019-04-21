from shapely.geometry import Point, Polygon
from shapely.affinity import scale, rotate
import math
try:
    long
except NameError:
    long = int

class Vector2:
    def __init__(self, param1, param2):
        if(isinstance(param1, (int, long, float)) and isinstance(param2, (int, long, float))):
            self.x = float(param1)
            self.y = float(param2)
        elif(isinstance(param1, Point) and isinstance(param2, Point)):
            self.x = param2.x - param1.x
            self.y = param2.y - param1.y
        else:
            raise ValueError("Vector Input Parameter Errors. Must be Numeric or Point")

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Vector2 (%.2lf, %.2lf)" % (self.x, self.y)

    def __copy__(self):
        return self.__class__(self.x, self.y)

    copy = __copy__

    def __len__(self):
        return 2

    def __getitem__(self, i):
        if(i > 1):
            raise IndexError("Vector2 index is 0 or 1.")
        return (self.x, self.y)[i]

    def __setitem__(self, i, value):
        if(i > 1):
            raise IndexError("Vector2 index is 0 or 1.")
        v = [self.x, self.y]
        v[i] = value
        self.x, self.y = l

    def __eq__(self, other):
        assert isinstance(other, Vector2)
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not self.__eq__(other)

    def __nonzero__(self):
        return self.x != 0 or self.y != 0

    def __neg__(self):
        return Vector2(-self.x, -self.y)

    def __add__(self, other):
        assert isinstance(other, Vector2)
        return Vector2(self.x + other.x, self.y + other.y)

    def __radd__(self, other):
        '''
        Vector + Vector = Vector
        Point + Vector = Point, point move along vector
        '''
        assert type(other) in [Vector2, Point]
        if isinstance(other, Vector2):
            return Vector2(self.x + other.x, self.y + other.y)
        if isinstance(other, Point):
            return Point(other.x + self.x, other.y + self.y)
    
    def __sub__(self, other):
        assert isinstance(other, Vector2)
        return Vector2(self.x - other.x, self.y - other.y)

    def __rsub__(self, other):
        assert isinstance(other, Vector2)
        return Vector2(other.x - self.x, other.y - self.y)

    def __mul__(self, other):
        assert isinstance(other, (int, long, float))
        return Vector2(self.x * other, self.y * other)

    __rmul__ = __mul__

    def __imul__(self, other):
        assert isinstance(other, (int, long, float))
        self.x *= other
        self.y *= other
        return self

    def __div__(self, other):
        assert type(other) in (int, long, float)
        return Vector2(self.x / other,self.y / other)

    def __abs__(self):
        return math.sqrt(self.x ** 2 + self.y ** 2)

    magnitude = __abs__

    def magnitude2(self):
        return self.x ** 2 + self.y ** 2

    def normalize(self):
        d = self.magnitude()
        if d:
            self.x /= d
            self.y /= d
        return self

    def normalized(self):
        d = self.magnitude()
        if d:
            return Vector2(self.x / d, self.y / d)
        return self.copy()

    def dot(self, other):
        assert isinstance(other, Vector2)
        return self.x * other.x + self.y * other.y

    def cross_m(self, other):
        assert isinstance(other, Vector2)
        return self.x * other.y - self.y * other.x

    def reflect(self, normal):
        assert isinstance(normal, Vector2)
        d = 2 * (self.x * normal.x + self.y * normal.y)
        return Vector2(self.x - d * normal.x, self.y - d * normal.y)

    def project(self, other):
        assert isinstance(other, Vector2)
        n = other.normalized()
        return self.dot(n)*n

    def angle(self, other):
        assert isinstance(other, Vector2)
        return math.acos(self.dot(other) / (self.magnitude()*other.magnitude()))

class LineSegment:
    def __init__(self, fpt, tpt):
        self.fpt = fpt
        self.tpt = tpt

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "LineSegment ( %s, %s )" % (str(self.fpt), str(self.tpt))

    def coords(self):
        return [self.fpt.coords[0], self.tpt.coords[0]]

    def length(self):
        return math.hypot(self.fpt.x - self.tpt.x, self.fpt.y - self.tpt.y)

    def height(self, pt):
        v2 = [pt.x - self.fpt.x, pt.y - self.fpt.y]
        v1 = [self.tpt.x - self.fpt.x, self.tpt.y  - self.fpt.y]
        cross = v1[0] * v2[1] - v1[1] * v2[0]
        height = math.fabs(cross) / math.hypot(v1[0], v1[1])
        return height

    def distance(self, pt):
        op = Point(self.fpt.x, self.fpt.y)
        dx, dy = self.tpt.x - op.x, self.tpt.y - op.y

        ratio = ((pt.x - op.x) * dx + (pt.y - op.y) * dy) / (dx * dx + dy * dy)
        if ratio > 1:
            op = self.tpt
        elif ratio > 0:
            op.x += ratio * dx
            op.y += ratio * dy

        dx = pt.x - op.x
        dy = pt.y - op.y

        return math.hypot(dx, dy)

    def distance2(self, pt):
        op = Point(self.fpt.x, self.fpt.y)
        dx, dy = self.tpt.x - op.x, self.tpt.y - op.y
        #print self
        ratio = ((pt.x - op.x) * dx + (pt.y - op.y) * dy) / (dx * dx + dy * dy)
        if ratio > 1:
            op = self.tpt
        elif ratio > 0:
            op.x += ratio * dx
            op.y += ratio * dy

        dx = pt.x - op.x
        dy = pt.y - op.y

        return dx * dx + dy * dy

    def vector(self):
        return Vector2(self.tpt.x - self.fpt.x, self.tpt.y - self.fpt.y)

class Circle:
    def __init__(self, center, radius):
        assert isinstance(center, Point)
        assert isinstance(radius, (int, long, float))
        self.center = center
        self.radius = radius

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Circle ( Center: %s, Radius: %s )" % (str(self.center), str(self.radius))

    @property
    def perimeter(self):
        return 2 * math.pi * self.radius

    @property
    def area(self):
        return math.pi * self.radius * self.radius

    def incircle(self, pt):
        return self.center.distance(pt) < self.radius
    
    @property
    def to_poly(self):
        return self.center.buffer(self.radius)

class Ellipse:
    def __init__(self, center, radii, rotation):
        self.center = center
        self.radii = radii
        self.rotation = rotation
        if(radii[0] < radii[1]):
            self.radii[0], self.radii[1] = radii[1], radii[0]
            self.rotation = [rotation[1], rotation[0]]

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Ellipse ( Center: %s, Radii: %s, Rotation %s)" % (str(self.center), str(self.radii), str(self.rotation))

    @property
    def e(self):
        a2 = self.radii[0] * self.radii[0]
        b2 = self.radii[1] * self.radii[1]
        return math.sqrt(1 - b2 / a2)

    @property
    def area(self):
        return math.pi * self.radii[0] * self.radii[1]

    @property
    def perimeter(self):
        '''
        Ramanujan
        '''
        a = self.radii[0]
        b = self.radii[1]
        h = math.pow((a - b) / (a + b), 2)
        p = math.pi * (a +b) * (1 + 3 * h / (10 + math.sqrt(4 - 3 * h)))
        return p

    @property
    def to_poly(self):
        cir = self.center.buffer(self.radii[0])
        e = self.radii[1] / self.radii[0]
        theta = math.atan2(self.rotation[0][1], self.rotation[0][0])
        return rotate(scale(cir, 1, e), theta, use_radians=True)


class Rectangle:
    def __init__(self, center, axes, extents):
        self.center = center
        self.axes = axes
        self.extents = extents
        if(extents[0] < extents[1]):
            self.extents = [extents[1], extents[0]]
            self.axes = [axes[1], axes[0]]

    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Rectangle ( Center: %s, Axes: %s, Extents: %s )" % (str(self.center), str(self.axes), str(self.extents))

    def __getitem__(self, i):
        if(i > 4 or i < 0):
            raise IndexError("Vector2 index is 0 ,1, 2, 3.")
        return self.coords[i]

    @property
    def perimeter(self):
        return (self.extents[0] + self.extents[1]) * 2

    @property
    def area(self):
        return self.extents[0] * self.extents[1]

    @property
    def length(self):
        return self.extents[0]

    @property
    def width(self):
        return self.extents[1]
    
    @property
    def to_poly(self):
        p1 = self.center + self.axes[0] * (0.5 * self.extents[0]) + self.axes[1] * (0.5 * self.extents[1])
        p2 = self.center + self.axes[0] * (-0.5 * self.extents[0]) + self.axes[1] * (0.5 * self.extents[1])
        p3 = self.center + self.axes[0] * (-0.5 * self.extents[0]) + self.axes[1] * (-0.5 * self.extents[1])
        p4 = self.center + self.axes[0] * (0.5 * self.extents[0]) + self.axes[1] * (-0.5 * self.extents[1])
        return Polygon([[p1.x, p1.y], [p2.x, p2.y], [p3.x, p3.y], [p4.x, p4.y]])
