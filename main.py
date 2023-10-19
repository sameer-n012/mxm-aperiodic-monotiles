import numpy as np
import itertools
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch
from HeeschSat.HexGrid import HexGrid

print('hi')

shape = np.array([
    [0, 0],
    [0, 1],
    [1, 0]
])

hg = HexGrid((6, 6))
ph = PolyhexagonHeesch(shape, 3, (6, 6))
ph.generate_variables()