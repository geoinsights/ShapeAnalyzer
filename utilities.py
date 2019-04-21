from shapely.geometry.base import BaseGeometry
from shapely.geometry import *
import pyvoronoi 

def dump_coords(geom):
    """Dump coordinates of a geometry into points set regardless of order"""
    points = []
    if not isinstance(geom, BaseGeometry):
        raise ValueError('Must be instance of a geometry class; found ' +
                         geom.__class__.__name__)
    elif geom.type in ('Point', 'LineString', 'LinearRing'):
        points = geom.coords[:]
        return points
    elif geom.type == 'Polygon':
        points.extend(geom.exterior.coords[:-1])
        for i in geom.interiors:
            points.extend(i.coords[:-1])
        return points
    elif geom.type.startswith('Multi') or geom.type == 'GeometryCollection':
        # Recursive call
        for part in geom:
            points.extend(dump_coords(part))
        return points
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))

def MaximumInscribedCircle(geom):
    '''
        Computer maximum inscribed circle of a polygon by voronoi of its boundary segments with pyvoronoi
    '''
    
    if not isinstance(geom, Polygon):
        raise ValueError('Must be instance of a Polygon class; found ' + geom.__class__.__name__)
    
    pv = pyvoronoi.Pyvoronoi(1)
    exterior = geom.exterior
    ptnum = len(exterior.coords)
    coords = exterior.coords[:]
    for i in range( ptnum - 1):
        j = (i+1) % ptnum
        pv.AddSegment([coords[i], coords[j]])

    for interior in geom.interiors:
        ptnum = len(interior.coords)
        coords = interior.coords[:]
        for i in range( ptnum - 1):
            j = (i+1) % ptnum
            pv.AddSegment([coords[i], coords[j]])

    pv.Construct()
    vertices = pv.GetVertices()

    radius = 0
    center = None
    for v in vertices:
        pt = Point(v.X, v.Y)
        if(geom.contains(pt)):
            r = pt.distance(geom.boundary)
            if(r > radius):
                radius = r
                center = pt
    return (center.x, center.y, radius)