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
        realRoots = r.real[abs(r.imag)<1e-5]
        for i in realRoots:
            if(0.01 < i and i < 0.99):
                return i, np.polyval(p1, i)

VLE_1atm = [7.44200335, -18.71029508, 18.98056828, -9.66489131, 2.95151852, 0]
Azeotrope = [1,0]
x, y = NumericalMethods.getPolynomialIntercept(VLE_1atm, Azeotrope)
print(str(x) + ', ' + str(y))

fig, ax = plt.subplots(1,1)
ax.plot([np.polyval(VLE_1atm, i) for i in np.linspace(0,1)], np.linspace(0,1))
ax.plot([np.polyval(Azeotrope, i) for i in np.linspace(0,1)], np.linspace(0,1))
ax.plot([0,x], [y,y], linestyle='dashed', c='black')
ax.plot([x,x], [0,y], linestyle='dashed', c='black')

plt.show()