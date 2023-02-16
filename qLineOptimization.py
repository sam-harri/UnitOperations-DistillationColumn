from typing import List
import numpy as np
from NumericalMethods.NumericalMetods import NumericalMethods
from math import sqrt

class qLineOptimazation:
    
    def distance(intersection : List[float], VLE : List[float]) -> float:
        x, y = intersection
        ydist = y - np.polyval(VLE, x)
        
        xt, yt = NumericalMethods.getPolynomialIntercept([y], VLE)
        xdist = xt - x
        return sqrt(xdist**2 + ydist**2)
    
    #unfinished, needs to get intersection point for given reflux ratio and xd
    def bestDistance(xd, R, VLE) -> float:
        qVals = np.linspace(0,1.5, num=1000)
        maxDistance = 0
        intersection = 0
        for q in qVals:
            distance = qLineOptimazation.distance(intersection, VLE)
            maxDistance = max(maxDistance, distance)
        return maxDistance