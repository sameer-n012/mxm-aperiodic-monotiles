import numpy as np
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch

if __name__ == '__main__':

    print('Test Starting...')

    t1 = np.array([
        [0, 0],
        [0, 1],
        [-1, 0]
    ])

    t2 = np.array([
        [0, 0],
        [0, 1]
    ])

    t3 = np.array([
        [0, -2],
        [0, -1],
        [0, 0],
        [0, 1],
        [-1, -2],
        [-1, 2]
    ])

    ph = PolyhexagonHeesch(t1, 2, (12, 12))
    ph.generate_variables()
    ph.construct_sat()
    ph.solve_sat()
    ph.write(directory='tests/out')

    print(len(ph.transforms))
    count = [0]*5
    for key in ph.transforms.keys():
        count[key[0]] += 1
    print('counts: ', count)

    exit()

    for k in range(0, 4):
        for s in range(max(6 * k, 6), 8 * k):
            ph = PolyhexagonHeesch(t1, k, (s, s))
            ph.generate_variables()
            ph.construct_sat()
            ph.solve_sat()
            ph.write(directory='tests/out')

            ph = PolyhexagonHeesch(t2, k, (s, s))
            ph.generate_variables()
            ph.construct_sat()
            ph.solve_sat()
            ph.write(directory='tests/out')

    print('Test Complete')
