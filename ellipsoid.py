import numpy as np
from numpy import linalg
from scipy import optimize
from random import random
from geom import Ellipse, Point, Vector2, Circle

class EnclosingEllipsoid:
    '''
    reference:
    https://github.com/minillinim/ellipsoid
    '''
    def getMinVolEllipse(self, P=None, tolerance=0.0001):
        """ Find the minimum volume ellipsoid which holds all the points
        
        Based on work by Nima Moshtagh
        http://www.mathworks.com/matlabcentral/fileexchange/9542
        and also by looking at:
        http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
        Which is based on the first reference anyway!
        
        Here, P is a numpy array of N dimensional points like this:
        P = [[x,y,z,...], <-- one point per line
             [x,y,z,...],
             [x,y,z,...]]
        
        Returns:
        (center, radii, rotation)
        
        """
        (N, d) = np.shape(P)
        d = float(d)
    
        # Q will be our working array
        Q = np.vstack([np.copy(P.T), np.ones(N)]) 
        QT = Q.T
        
        # initializations
        err = 1.0 + tolerance
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > tolerance:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            M = np.diag(np.dot(QT , np.dot(linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipse 
        center = np.dot(P.T, u)
    
        # the A matrix for the ellipse
        A = linalg.inv(
                       np.dot(P.T, np.dot(np.diag(u), P)) - 
                       np.array([[a * b for b in center] for a in center])
                       ) / d
                       
        # Get the values we'd like to return
        U, s, rotation = linalg.svd(A)
        radii = 1.0/np.sqrt(s)
        V1 = Vector2(rotation[0][0], rotation[0][1])
        V2 = Vector2(rotation[1][0], rotation[1][1])
        return Ellipse(Point(center), radii, [V1, V2])

class LSqEllipse:
    '''
    reference : https://github.com/bdhammel/least-squares-ellipse-fitting
    '''
    def fit(self, data):
        x, y = np.asarray(data, dtype=float)
        cx, cy = x.mean(), y.mean()
        x = x - cx
        y = y - cy
        #Quadratic part of design matrix [eqn. 15] from (*)
        D1 = np.mat(np.vstack([x**2, x*y, y**2])).T
        #Linear part of design matrix [eqn. 16] from (*)
        D2 = np.mat(np.vstack([x, y, np.ones(len(x))])).T
        
        #forming scatter matrix [eqn. 17] from (*)
        S1 = D1.T*D1
        S2 = D1.T*D2
        S3 = D2.T*D2  
        
        #Constraint matrix [eqn. 18]
        C1 = np.mat('0. 0. 2.; 0. -1. 0.; 2. 0. 0.')

        #Reduced scatter matrix [eqn. 29]
        M=C1.I*(S1-S2*S3.I*S2.T)

        #M*|a b c >=l|a b c >. Find eigenvalues and eigenvectors from this equation [eqn. 28]
        eval, evec = np.linalg.eig(M) 

        # eigenvector must meet constraint 4ac - b^2 to be valid.
        cond = 4*np.multiply(evec[0, :], evec[2, :]) - np.power(evec[1, :], 2)
        a1 = evec[:, np.nonzero(cond.A > 0)[1]]
        
        #|d f g> = -S3^(-1)*S2^(T)*|a b c> [eqn. 24]
        a2 = -S3.I*S2.T*a1
        
        # eigenvectors |a b c d f g> 
        self.coef = np.vstack([a1, a2])
       
        #eigenvectors are the coefficients of an ellipse in general form
        #a*x^2 + 2*b*x*y + c*y^2 + 2*d*x + 2*f*y + g = 0 [eqn. 15) from (**) or (***)
        a = self.coef[0,0]
        b = self.coef[1,0]/2.
        c = self.coef[2,0]
        d = self.coef[3,0]/2.
        f = self.coef[4,0]/2.
        g = self.coef[5,0]
        
        #finding center of ellipse [eqn.19 and 20] from (**)
        x0 = (c*d-b*f)/(b**2.-a*c)
        y0 = (a*f-b*d)/(b**2.-a*c)
        
        #Find the semi-axes lengths [eqn. 21 and 22] from (**)
        numerator = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        denominator1 = (b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        denominator2 = (b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        width = np.sqrt(numerator/denominator1)
        height = np.sqrt(numerator/denominator2)

        # angle of counterclockwise rotation of major-axis of ellipse to x-axis [eqn. 23] from (**)
        # or [eqn. 26] from (***).
        phi = .5*np.arctan((2.*b)/(a-c))
        V1 = Vector2(1, np.tan(phi)).normalize()
        V2 = Vector2(-V1.y, V1.x)
        x0 += cx
        y0 += cy
        return Ellipse(Point(x0, y0),[width, height],[V1, V2])

class LSqCircle:
    '''
    reference : 
    [1]https://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
    [2]https://gist.github.com/lorenzoriano/6799568
    '''
    def fit(self, data):
        x_m = np.mean(data[0])
        y_m = np.mean(data[1])
        center_estimate = x_m, y_m
        center, ier = optimize.leastsq(self._f, center_estimate, args=(data[0],data[1]))
        xc, yc = center
        Ri       = self._calc_R(data[0], data[1], *center)
        R        = Ri.mean()
        residu   = np.sum((Ri - R)**2)
        return Circle(Point(xc, yc),R)
    
    def _calc_R(self, x,y, xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((x-xc)**2 + (y-yc)**2)

    def _f(self, c, x, y):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = self._calc_R(x, y, *c)
        return Ri - Ri.mean()
    
ECELL = EnclosingEllipsoid()
LSELL = LSqEllipse()
LSCIR = LSqCircle()
