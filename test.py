from ReactorConstants import ReactorConstants
import numpy as np
import matplotlib.pyplot as plt
from typing import List


VLE_1atm_test : List[float] = ReactorConstants.VLE_1atm_test
AzeotropicPoly : List[float] = [1,0]

x = np.linspace(0,1)
y = [np.polyval(VLE_1atm_test, i) for i in x]

x2 = np.linspace(0,1,2)
y2 = [np.polyval(AzeotropicPoly, i) for i in x2]

fig, ax = plt.subplots(1,1)
ax.plot(x,y)
ax.plot(x2,y2)
plt.show()

