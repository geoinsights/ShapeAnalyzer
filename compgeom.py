import math
from shapely.geometry import Point
from geom import *

class GeometryUtility():
    def ccw(self, p1, p2, p3):
        return (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x)

    def area2(self, p1, p2, p3):
        return math.fabs(self.ccw(p1, p2, p3))

    def graham_scan(self,points):
        num = len(points)
        cpt = Point(0.0,0.0)
        for i in range(num-2):
            if self.ccw(points[i], points[i+1], points[i+2]) != 0:
                cptx = (points[i].x + points[i+1].x + points[i+2].x) / 3
                cpty = (points[i].y + points[i+1].y + points[i+2].y) / 3
                cpt = Point(cptx, cpty)
                break
        angles = []
        distances = []
        for i in range(num):
            distances.append(cpt.distance(points[i]))
            dx = points[i].x - cpt.x
            dy = points[i].y - cpt.y
            th = math.atan2(dy, dx)
            angles.append(th)
        idx = range(num)
        data = zip(idx, points, angles, distances)
        data = sorted(data, key = lambda it: it[2])
        i = 1
        convex_runs = 0
        kk = 0
        while (convex_runs < num):
            kk +=1
            mod = len(data)
            curr = data[i][1]
            forw = data[(i+1) % mod][1]
            back = data[(i-1) % mod][1]
            if self.ccw(back, curr, forw) <= 0:
                data.pop(i)
                convex_runs = 0
            else:
                convex_runs += 1
            mod = len(data)
            i = (i+1) % mod
        return list(zip(*data))[1]

GU = GeometryUtility()

class RotatingCaliper:
    def __init__(self, points):
        self._points = points
        self._N = len(points)
        self._antipodals = None
        self._calipers()

    def _next(self, i):
        return (i+1) % self._N

    def _calipers(self):
        q = 1
        antipodals = {}
        N = len(self._points)
        for i in range(N):
            antipodals[i] = [-1,-1]
        next = self._next
        PSet = self._points
        for i in range(N):
            j = next(i)
            while GU.area2(PSet[i], PSet[j], PSet[next(q)]) > GU.area2(PSet[i], PSet[j], PSet[q]):
                q = next(q)
            antipodals[j][0] = q
            antipodals[i][1] = q
        self._antipodals = antipodals

    def AntiPodals(self):
        antis = []
        PSet = self._points
        for k, v in self._antipodals.items():
            for i in range(v[0], v[1]+1):
                if [PSet[k], PSet[i]] in antis or [PSet[i], PSet[k]] in antis:
                    continue
                antis.append([PSet[k], PSet[i]])
        return antis

    def Diameter(self):
        antis = self.AntiPodals()
        d = -1
        support = []
        for ppair in antis:
            td = ppair[0].distance(ppair[1])
            if td > d:
                d = td
                support = ppair
        return d, support

    def Width(self):
        w = float("inf")
        support = []
        for i in range(self._N):
            fpt = self._points[i]
            tpt = self._points[(i+1) % self._N]
            lseg = LineSegment(fpt, tpt)
            j = self._antipodals[i][-1]
            pt = self._points[j]
            td = lseg.height(pt)
            if td < w:
                w = td
                support = [fpt, tpt, pt]
        return w, support