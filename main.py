import numpy as np
import itertools
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch
from HeeschSat.HexGrid import HexGrid

print('hi')

# shape = np.array([
#     [0, 0],
#     [0, 1],
#     [1, 0]
# ])

# hg = HexGrid((6, 6))
# ph = PolyhexagonHeesch(shape, 3, (6, 6))
# ph.generate_variables()

six_hex = np.array([
    [0,0],
    [0,1],
    [0,-1],
    [-1,2],
    [0,-2],
    [-1,-2]
])

ph = PolyhexagonHeesch(six_hex, 1, (10, 10))
ph.generate_variables()
print('generated')
ph.construct_sat()
print(ph.sat)
ph.solve_sat()

m = ph.model

minimum = 10000
for i in m:
    if i <= 0:
        continue
    print(i)
    if i < minimum:
        minimum = i

for k, v in ph.transforms.items():
    if v[0] == minimum:
        print(k)


print('a')
print(m)


