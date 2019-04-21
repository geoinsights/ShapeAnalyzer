from shapely.affinity import rotate, translate, scale
from shapely.geometry import Polygon, Point
import numpy as np

from mincircle import make_circle
from utilities import *
from geom import *
from compgeom import *
from ellipsoid import LSCIR, LSELL, ECELL

class Shape(Polygon, MultiPolygon):
    '''
    shape by polygon  
    '''
    @property
    def coords(self):
        '''
        shape boundary coordinates list without order
        '''
        return dump_coords(self)

    @property
    def boundary_coords(self):
        return self.boundary.coords[:-1]

    @property
    def to_points(self):
        pts = []
        for coord in self.coords:
            pts.append(Point(coord[0], coord[1]))
        return pts

    @property
    def distance_angle_matrix(self):
        coords = np.array(self.coords)
        dmat = [np.hypot((c - coords)[:,0],(c - coords)[:,1]) for c in coords]
        dmat = np.array(dmat)
        amat = [np.arctan2((coords - c)[:,1], (coords - c)[:,0]) for c in coords]
        amat = np.array(amat)

        amat[amat <0] += np.pi *2 
        amat = amat / np.pi * 180
        return dmat, amat
        
    '''
    shape transform: 
        translate, rotate, scale.
    '''
    def translate(self, xoff, yoff):
        return Shape(translate(self, xoff, yoff))

    def rotate(self, theta):
        return Shape(rotate(self, theta, 'centroid'))

    def scale(self, fact):
        return Shape(scale(self, xfact=fact,yfact=fact,origin='centroid'))

    '''
    shape bounding geometry: 
        bbox, convex hull, 
        minimum area rectangle, minimum width rectangle, 
        minimum enclosing circle.
        minimum enclosing ellipse
    '''
    @property
    def bbox(self):
        rect = self.envelope
        coords = rect.exterior.coords
        center = rect.centroid
        V1 = Vector2(Point(coords[0]), Point(coords[1])).normalize()
        V2 = Vector2(Point(coords[1]), Point(coords[2])).normalize()
        ext1 = Point(coords[0]).distance(Point(coords[1]))
        ext2 = Point(coords[1]).distance(Point(coords[2]))
        return Rectangle(center, [V1, V2], [ext1, ext2])

    @property
    def chull(self):
        return Shape(self.convex_hull)
    
    @property
    def minimum_area_rect(self):
        rect = self.minimum_rotated_rectangle
        coords = rect.exterior.coords
        center = rect.centroid
        V1 = Vector2(Point(coords[0]), Point(coords[1])).normalize()
        V2 = Vector2(Point(coords[1]), Point(coords[2])).normalize()
        ext1 = Point(coords[0]).distance(Point(coords[1]))
        ext2 = Point(coords[1]).distance(Point(coords[2]))
        return Rectangle(center, [V1, V2], [ext1, ext2])

    @property
    def minimum_width_rect(self):
        hull = GU.graham_scan(self.to_points)
        rc = RotatingCaliper(hull)
        width, support = rc.Width()
        min0, max0 = 0, 0
        V0 = Vector2(support[1].x - support[0].x, support[1].y - support[0].y)
        V0.normalize()
        Vh = Vector2(-V0.y, V0.x)
        Vh.normalize()
        for i in range(len(hull)):
            Vi = Vector2(hull[i].x - support[0].x, hull[i].y - support[0].y)
            dot = V0.dot(Vi)
            if dot < min0:
                min0 = dot
            elif dot > max0:
                max0 = dot
        #area = (max0 - min0) * width
        extx = max0 - min0
        exty = width
        cx = support[0].x + (min0 + max0) / 2.0 * V0.x + width / 2.0 * Vh.x
        cy = support[0].y + (min0 + max0) / 2.0 * V0.y + width / 2.0 * Vh.y
        center = Point(cx, cy)
        rect = Rectangle(center, [V0, Vh], [extx, exty])
        return rect

    @property
    def minimum_enclosing_circle(self):
        cx, cy, R = make_circle(self.coords)
        return Circle(Point(cx, cy), R)

    @property
    def minimum_enclosing_ellipse(self):
        coords = self.coords
        P = np.reshape(coords, (len(coords), 2))
        return ECELL.getMinVolEllipse(P)

    '''
    fitting shape coords with model shape: circle, ellipse
    '''
    @property
    def fitting_ellipse(self):
        coords = np.array(self.coords)
        return LSELL.fit([coords[:,0], coords[:,1]])
        
    @property
    def fitting_circle(self):
        coords = np.array(self.coords)
        return LSCIR.fit([coords[:,0], coords[:,1]])

    '''
    maximum inscribed circle
    '''
    @property
    def maximum_inscribed_circle(self):
        cx, cy, R = MaximumInscribedCircle(self)
        return Circle(Point(cx, cy), R)
    
    '''
    centers of shape
    '''
    @property
    def coords_mean(self):
        '''
        mean coordinates: [mean(X), mean(Y)]
        '''
        return Point(np.mean(self.coords,0))
    
    @property
    def area_centroid(self):
        '''
        centroid of area: weighted sum of centroids of triangles in the area triangulation
        '''
        return self.centroid
    
    @property
    def boundary_centroid(self):
        '''
        centroid of boundary: weighted sum of centroids of each segment of boundary
        '''
        return self.boundary.centroid
    
    @property
    def shape_eigen(self):
        '''
        coords eigen value and eigne vector, for orientation of shape
        '''
        coords = np.asarray(self.coords) # remove the repeated first point
        center = self.area_centroid
        coords -= [center.x, center.y]
        A = np.cov(coords.T)
        eigval, eigvet = np.linalg.eig(A)
        order = eigval.argsort()[::-1]
        eigval = eigval[order]
        eigvet = eigvet[:,order]
        V1 = Vector2(*eigvet[:,0])
        V2 = Vector2(*eigvet[:,1])
        return eigval, (V1, V2)

    @property
    def least_inertia(self):
        '''
        orientation of shape, angle with x-axis +
        '''
        eigval, eigvet = self.shape_eigen
        theta = np.arctan2(eigvet[0].y, eigvet[0].x)
        if (theta<0):
            theta += np.pi
        return theta
    
    @property
    def principal_axes_ratio(self):
        eigval, eigvet = self.shape_eigen
        return np.sqrt(eigval[1] / eigval[0])

    def confidence_ellipse(self, fact=1.0):
        eigval, eigvet = self.shape_eigen
        center = self.area_centroid
        return Ellipse(center, np.sqrt(eigval) * fact, eigvet)

    @property
    def ellipse_var(self):
        coords = np.asarray(self.coords)
        center = self.area_centroid
        coords -= [center.x, center.y] 
        N = len(coords)
        a = np.sum(coords[:,0] * coords[:,0]) / N
        b = 2 * np.sum(coords[:,0] * coords[:,1]) / N
        c = np.sum(coords[:,1] * coords[:,1]) / N
        cor = np.asarray([[a,b/2],[b/2,c]])
        
        dis = [c.dot(np.linalg.inv(cor)).dot(c.T) for c in coords]
        dis = np.sqrt(dis)
        mu_r = np.mean(dis)
        std_r = np.std(dis)
        return std_r / mu_r
    
    @property
    def elongation(self):
        rect = self.minimum_width_rect
        elg = 1 - rect.width / rect.length
        return elg
        
    @property
    def rectangularity(self):
        area = self.area
        area_r = self.minimum_area_rect.area
        return area / area_r
    
    @property
    def convexity(self):
        perimeter = self.boundary.length
        perimeter_ch = self.convex_hull.boundary.length
        return perimeter_ch / perimeter
    
    @property
    def solidity(self):
        area = self.area
        area_ch = self.convex_hull.area
        return area / area_ch
    
    @property
    def circularity_area(self):
        area = self.area
        perimeter = self.boundary.length
        return 4 * np.pi * area / (perimeter **2)
    
    @property
    def circularity_radial(self):
        center = self.area_centroid
        coords = self.coords
        radials = [center.distance(Point(c[0], c[1])) for c in coords]
        mu_r = np.mean(radials)
        std_r = np.std(radials)
        return std_r / mu_r
    
    @property
    def euler_number(self):
        return 1 - len(self.interiors)
    
    @property
    def holes_ratio(self):
        area = self.area
        area_h = 0.0
        for interior in self.interiors:
            area_h += Polygon(interior).area
        return area_h / area
    
    @property
    def complex_coords(self):
        center = self.area_centroid
        coords = np.asarray(self.coords)
        coords -= [center.x, center.y]
        cp = [np.complex(c[0], c[1]) for c in coords]
        return cp

    @property
    def centroid_distance(self):
        center = self.area_centroid
        points = self.to_points
        dis = [center.distance(p) for p in points]
        return dis

    @property
    def tangent_function(self):
        '''
        step function of line segments, turn angle function
        '''
        coords = dump_coords(self.boundary)[:-1]
        length = self.boundary.length
        N = len(coords)
        L = 0
        result = []
        for i in range(N):
            j = (i+1) % N
            pi = Point(coords[i])
            pj = Point(coords[j])
            ti = L / length
            L += pi.distance(pj)
            tj = L / length
            theta = np.angle(np.complex(pj.x-pi.x,pj.y-pi.y))
            if(theta < 0):
                theta += 2 * np.pi
            result.append([ti, theta])
            result.append([tj, theta])
        return result
    
    @property
    def area_function(self):
        coords = self.boundary_coords
        N = len(coords)
        centroid = self.area_centroid
        result = []
        for i in range(N):
            j = (i+1) % N
            pi = coords[i]
            pj = coords[j]
            tri = Polygon([centroid.coords[0], pi, pj])
            result.append(tri.area)
        return result
    
    def tar_signature(self, t):
        coords = self.exterior.coords[:-1]
        if( self.exterior.is_ccw is False):
            coords.reverse()
        N = len(coords)
        if(t > round(N / 2)):
            raise ValueError("parameter t cannot be larger than N / 2")
        results = []
        for j in range(N):
            i = j - t
            if(i < 0):
                i += N
            k = (j + t) % N
            xi,yi = coords[i]
            xj,yj = coords[j]
            xk,yk = coords[k]
            sarea = (xi*yj - xi*yk - xj*yi + xk*yi + xj*yk - xk*yj) * 0.5
            results.append(sarea)
        return results

    def discrete_curve_evolution(self):
        pass
    
    def shape_dft(self, N):
        pass

    def grid_matrix(self):
        pass

    def polar_matrix(self):
        pass

    def chain_code(self, n):
        pass

    def beam_angle(self, n):
        pass
    



