import numpy as np
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch
from HeeschSat.PolykiteHeesch import PolykiteHeesch
if __name__ == '__main__':

    print('Test Starting...')

    t0 = np.array([
        [0, 0],
    ])

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

    t4 = np.array([
        [0, 0],
        [0, -1],
        [0, -2],
        [-1, 1],
        [1, -1],
        [2, -1],
        [2, 0]
    ])

    tk1 = np.array([
        [0, 0, 0],
        [1, 0, 2]
    ])

    ph = PolykiteHeesch(tk1, coronas=0)
    ph.generate_variables()
    ph.construct_sat()
    ph.solve_sat()
    ph.write(directory='tests/out', plot=False)

    # Count transforms in each corona
    print(len(ph.transforms))
    count = [0]*(ph.k_cor+1)
    for key in ph.transforms.keys():
        count[key[0]] += 1
        # if key[0] == 1 and (key[1] > 7 or key[2] > 7):
        #     print(key)
        #     print(ph.transforms[key])
        #     break
    print('counts: ', count)
    print('time: ', ph.times[-1] - ph.times[0])

    print('Test Complete')

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
