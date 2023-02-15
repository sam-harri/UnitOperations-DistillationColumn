from typing import List
import numpy as np
import matplotlib.pyplot as plt

class NumericalMethods:
    
    @staticmethod
    def getPolynomialIntercept(p1 : List[float], p2 : List[float]):
        
        if(len(p2)<len(p1)):
            poly = p1
            add = p2
        else:
            poly = p2
            add = p1
            
        add, poly = add[::-1], poly[::-1]
        for i in range(len(add)):
            poly[i] -= add[i]
        
        poly = poly[::-1]
        
        r = np.roots(poly)
        realRoots = r.real[abs(r.imag)<1e-4]
        for i in realRoots:
            if(0.01 < i and i < 0.99):
                return i, np.polyval(p1, i)
    
    @staticmethod
    def staircaseDataAzeotrope(p1 : List[float], p2 : List[float], startPoint : List[float], xd : float):
        #p1 = VLE
        #p2 = Azeotrope Line
        data : List[List[float]] = []
        x, y = startPoint
        
        while(x < xd):
            data.append([x,y])
            y = np.polyval(p1,x)
            data.append([x,y])
            x,y = NumericalMethods.getPolynomialIntercept([y],p2)
        data.append([x,y])

        return data
    
    @staticmethod
    def staircaseData1atm(oLine : List[float], sLine : List[float], VLELine : List[float], xd: float, xb : float):
        xi, yi = NumericalMethods.getPolynomialIntercept(sLine, oLine)
        data : List[List[float]]= []
        x, y = xb, xb
        
        while(x < xd):
            data.append([x,y])
            y = np.polyval(VLELine,x)
            data.append([x,y])
            if(y>=yi): #in operating section now
                x,y = NumericalMethods.getPolynomialIntercept([y],oLine)
            else: #in stripping
                x,y = NumericalMethods.getPolynomialIntercept([y], sLine)
        
        if(y>xd+0.01):
            data = data[:-1]
        else:
            data.append([x,y])
        return data
    
    @staticmethod
    def staircaseData10atm(oLine : List[float], sLine : List[float], VLELine : List[float], xd: float, xb : float):
        xi ,yi = NumericalMethods.getPolynomialIntercept(oLine, sLine)
        data : List[List[float]] = []
        x, y = xb, xb
        
        while(x > xd):
            data.append([x,y])
            y = np.polyval(VLELine,x)
            data.append([x,y])
            if(y<=yi): #in operating section now
                x,y = NumericalMethods.getPolynomialIntercept([y],oLine)
            else: #in stripping
                x,y = NumericalMethods.getPolynomialIntercept([y], sLine)
        
        if(y+0.01<xd):
            data = data[:-1]
        else:
            data.append([x,y])
        return data
    
    def qLinePolynomial(q, zf):
        if q == 1:
            q=0.999
        return [(q)/(q-1),(-zf)/(q-1)]
    
    def operatingLinePolynomial(refluxRatio : float, xd : float):
        return [(refluxRatio)/(refluxRatio+1),(xd)/(refluxRatio+1)]
    
    def strippingLinePolynomial(qLine : List[float], oLine : List[float], xb : float):
        x1, y1 = xb, xb
        x2, y2 = NumericalMethods.getPolynomialIntercept(qLine, oLine)
        m = (y2-y1)/(x2-x1)
        b = xb - (m*xb)
        return [m,b]