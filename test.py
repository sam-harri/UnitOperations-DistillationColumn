from ReactorConstants import ReactorConstants
import numpy as np
import matplotlib.pyplot as plt
from typing import List


tmp = [[1,1], [0.456, 0.01]]
tmp2 = [[130], [32.5]]

print(np.linalg.solve(tmp, tmp2))
