from typing import List
import numpy as np
import matplotlib.pyplot as plt

class NumericalMethods:
    
    @staticmethod
    def getPolynomialIntercept(p1 : List[float], p2 : List[float]):# -> tuple(float,float):
        
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
        realRoots = r.real[abs(r.imag)<1e-3]
        for i in realRoots:
            if(0.01 < i and i < 0.99):
                return i, np.polyval(p1, i)
        raise Exception(f"No Valid Intersection Between Polynomials {p1} and {p2}")
    
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
    def staircaseData1atm(q, zf, R, VLELine : List[float], xd, xb, iterations):# -> tuple(int,List[List[float]],List[List[float]]):
        if(iterations>5):
            raise Exception("Non Convering 1atm Column Parameters")
        qLine = NumericalMethods.qLinePolynomial(q, zf)
        oLine = NumericalMethods.operatingLinePolynomial(R, xd)
        if(qLine == None or oLine == None):
            raise Exception("Invalide qLine or oLine")
        sLine = NumericalMethods.strippingLinePolynomial(qLine, oLine, xb)
        if(sLine == None):
            raise Exception("Invalid sLine based on qLine and oLine")
        xi, yi = NumericalMethods.getPolynomialIntercept(sLine, oLine)
        if(yi > np.polyval(VLELine, xi)):
            raise Exception("1atm Pinch Point Above VLE Line")
        data : List[List[float]]= []
        x, y = xb, xb
        steps = 0
        while(x < xd):
            data.append([x,y])
            y = np.polyval(VLELine,x)
            data.append([x,y])
            if(y>=yi): #in operating section now
                x,y = NumericalMethods.getPolynomialIntercept([y],oLine)
            else: #in stripping
                x,y = NumericalMethods.getPolynomialIntercept([y], sLine)
            steps +=1

        data.append([y,y])
        
        if(data[-1][1] > xd+0.02):
            # print("is " + str(data[-1][1]) + " but should be " + str(xd))
            return NumericalMethods.staircaseData1atm(q, zf, R, VLELine, 0.5 * (xd + data[-1][1]), xb, iterations + 1)
        elif(data[-1][1] < xd-0.02):
            # print("is " + str(data[-1][1]) + "but should be " + str(xd))
            return NumericalMethods.staircaseData1atm(q, zf, R, VLELine, 0.5 * (xd + data[-1][1]), xb, iterations+1)
        else:
            # print("is " + str(data[-1][1]) + " but should be " + str(xd) + " : Done")
            lines = [qLine, oLine, sLine]
            if(data[-1][0] > xd+0.01):
                data.pop()
                data.pop()
            return steps, data, lines
    @staticmethod
    def staircaseData10atm(q, zf, R, VLELine, xd, xb, iterations):
        if(iterations>5):
            raise Exception("Non Convering 10atm Column Parameters")
        qLine = NumericalMethods.qLinePolynomial(q, zf)
        oLine = NumericalMethods.operatingLinePolynomial(R, xd)
        if(qLine == None or oLine == None):
            raise Exception("Invalide qLine or oLine")
        sLine = NumericalMethods.strippingLinePolynomial(qLine, oLine, xb)
        if(sLine == None):
            raise Exception("Invalid sLine based on qLine and oLine")
        xi, yi = NumericalMethods.getPolynomialIntercept(sLine, oLine)
        if(yi < np.polyval(VLELine, yi)):
            raise Exception("10atm Pinch Point Below VLE Line")
        data : List[List[float]] = []
        x, y = xb, xb
        steps = 0
        while(x > xd):
            data.append([x,y])
            y = np.polyval(VLELine,x)
            data.append([x,y])
            if(y<=yi): #in operating section now
                x,y = NumericalMethods.getPolynomialIntercept([y],oLine)
            else: #in stripping
                x,y = NumericalMethods.getPolynomialIntercept([y], sLine)
            steps +=1

        data.append([y,y])
        
        if(data[-1][1] > xd+0.02):
            # print("is " + str(data[-1][1]) + " but should be " + str(xd))
            return NumericalMethods.staircaseData10atm(q, zf, R, VLELine, 0.5 * (xd + data[-1][1]), xb, iterations+1)
        elif(data[-1][1] < xd-0.02):
            # print("is " + str(data[-1][1]) + "but should be " + str(xd))
            return NumericalMethods.staircaseData10atm(q, zf, R, VLELine, 0.5 * (xd + data[-1][1]), xb, iterations+1)
        else:
            # print("is " + str(data[-1][1]) + " but should be " + str(xd) + " : Done")
            lines = [qLine, oLine, sLine]
            
            return steps, data, lines
    
    def qLinePolynomial(q, zf) -> List[float]:
        if q == 1:
            q=0.999
        return [(q)/(q-1),(-zf)/(q-1)]
    
    def operatingLinePolynomial(refluxRatio, xd) -> List[float]:
        return [(refluxRatio)/(refluxRatio+1),(xd)/(refluxRatio+1)]
    
    def strippingLinePolynomial(qLine : List[float], oLine : List[float], xb) -> List[float]:
        x1, y1 = xb, xb
        x2, y2 = NumericalMethods.getPolynomialIntercept(qLine, oLine)
        if(x2==None or y2 ==None):
            raise Exception("No Valid Intersection Bewteen the qLine and Operating Line")
        m = (y2-y1)/(x2-x1)
        b = xb - (m*xb)
        return [m,b]
    
    def isValidIntersection(intersection : List[float], VLE : List[float]):
        x,y = intersection
        if(y > np.polyval(VLE, x) or y < x):
            return False
        return True