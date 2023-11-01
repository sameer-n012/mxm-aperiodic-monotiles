import numpy as np
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch

if __name__ == '__main__':

    print('Test Starting...')

    tile = np.array([
        [0, 0],
        [0, 1],
        [-1, 0]
    ])

    ph = PolyhexagonHeesch(tile, 2, (12, 12))
    ph.generate_variables()
    ph.construct_sat()
    ph.solve_sat()
    ph.write(directory='tests/out')
    print('Test Complete')
