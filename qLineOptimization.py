from typing import List
import numpy as np
from NumericalMethods.NumericalMetods import NumericalMethods
from ReactorConstants import ReactorConstants
from math import sqrt
import matplotlib.pylab as plt

class qLineOptimazation:
    
    @staticmethod
    def distance(intersection : List[float], VLE : List[float]) -> float:
        x, y = intersection
        ydist = np.polyval(VLE, x)
        xdist = NumericalMethods.getPolynomialIntercept(VLE, [y])[1]
        return sqrt(ydist**2 + xdist**2)
    
    #implement qLineOptimization.distance and NumericalMethods.isValidIntersection methods
    @staticmethod
    def optimized_q1atm(zf : float, R : float, xd : float, VLE : List[float]) -> float:
        
        q_arr = np.linspace(-1,2,301)
        bestq = [0,0] #q, distance (random small distance instanced first)
        qdata = []
        distancedata = []

        
        for q in q_arr:
            qLine = NumericalMethods.qLinePolynomial(q, zf)
            oLine = NumericalMethods.operatingLinePolynomial(R, xd)
            intersection = NumericalMethods.getPolynomialIntercept(qLine, oLine)
            if(intersection == None):
                continue
            if(NumericalMethods.isValidIntersection(intersection, VLE)):
                distance = qLineOptimazation.distance(intersection, VLE)
                qdata.append(q)
                distancedata.append(distance)
                if(distance > bestq[1]):
                    bestq = [q, distance]
        return bestq, qdata, distancedata

bestq, qdata, distancedata = qLineOptimazation.optimized_q1atm(0.25,1.35,0.454,ReactorConstants.VLE_1atm)
fig, ax = plt.subplots(1,1)
ax.plot(qdata, distancedata)
plt.show()

