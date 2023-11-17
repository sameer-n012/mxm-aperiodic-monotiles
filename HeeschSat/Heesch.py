from abc import ABC, abstractmethod
import numpy as np
import itertools
from pysat.solvers import Solver
from time import time
import matplotlib.colors as mcolors
from multiprocessing import Process, Manager, Pool, cpu_count


def define_colors():
    lst = sorted(mcolors.CSS4_COLORS,
                 key=lambda c:
                 tuple(
                     mcolors.rgb_to_hsv(mcolors.to_rgb(c))
                 ))
    lst.remove('black')
    return lst, [0, 14, 33, 53, 64, 85, 122, 132]


class Heesch(ABC):
    # plot_colors = ['k', 'b', 'r', 'y', 'g', 'c', 'm']
    plot_colors, c_majors = define_colors()
    c_spacing = 3
    num_timestamps = 14

    def __init__(self):
        self.k_cor = None
        self.grid = None
        self.transforms = {}
        self.cells = {}
        self.rotation_matrices = []
        self.shape = None
        self.shape_size = None
        self.shape_rad = None
        self.sat = None
        self.model = None
        self.num_clauses = 0
        self.times = [0.0] * Heesch.num_timestamps

    def generate_variables(self):

        self.times[0] = time()

        mid = int(self.grid.size[0] / 2), int(self.grid.size[1] / 2)

        max_transforms = (self.k_cor + 1) * self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)

        # include single 0-corona transform as original shape
        translate_mat = np.empty(self.shape_size, dtype=int)
        translate_mat[:, 0] = mid[0]
        translate_mat[:, 1] = mid[1]
        transform = self.shape + translate_mat
        halo = np.unique([i for c in transform for i in self.grid.haloIdx(c)], axis=0)
        key = (0, mid[0], mid[1], 0)
        self.transforms[key] = self.get_transform_idx(key), transform, halo

        # get the used rotation matrices (some might result in rotational symmetry)
        rotate_indices = self.check_rotational_symmetry()

        # transform variables (if a transform is used)
        offset = self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)
        for idx, val in enumerate(itertools.product(range(1, self.k_cor + 1),
                                                    range(0, self.grid.size[0]),
                                                    range(0, self.grid.size[1]),
                                                    range(0, len(self.rotation_matrices))
                                                    )):
            if val[0] == 0:
                continue

            if rotate_indices[val[3]] == 0:
                continue

            translate_mat = np.empty(self.shape_size, dtype=int)
            translate_mat[:, 0] = val[1]
            translate_mat[:, 1] = val[2]
            rotate_mat = self.rotation_matrices[val[3]]
            transform = np.matmul(rotate_mat, self.shape.T).T + translate_mat

            if max(abs(val[1] - mid[0]), abs(val[2] - mid[1])) > (val[0] + 1) * (self.shape_rad + 0):
                continue
            if max(abs(val[1] - mid[0]), abs(val[2] - mid[1])) < (val[0] - 1):
                continue

            # only include transforms that fall within the bounds of the grid
            if self.grid.in_bounds(transform):
                halo = np.unique([i for c in transform for i in self.grid.haloIdx(c)], axis=0)
                self.transforms[tuple(val)] = idx + offset + 1, transform, halo

        self.times[1] = time()

        # cell variables (if a cell is taken)
        offset = max_transforms + 1
        for idx, i in enumerate(self.grid.indices()):
            self.cells[tuple(i)] = offset + idx

        self.times[2] = time()

    def construct_sat(self) -> Solver:
        """
        Note that a, b => c | d === -a | -b | c | d
        """

        s = Solver(name="g3")

        self.times[3] = time()

        # 0-corona always used
        k_0 = self.transforms[(0, int(self.grid.size[0] / 2), int(self.grid.size[1] / 2), 0)]
        s.add_clause([k_0[0]])
        self.num_clauses += 1

        self.times[4] = time()

        # if a transform is used, its cells are used
        for k, v in self.transforms.items():

            # ad each cell in the transform
            for c in v[1]:
                s.add_clause([-v[0], self.cells[tuple(c)]])
                self.num_clauses += 1

        self.times[5] = time()

        # if a cell is used, some transform must use it
        for k1, v1 in self.cells.items():
            lst = [-v1]

            # check if the cell is in the transform for
            # each transform
            for k2, v2 in self.transforms.items():
                if list(k1) in v2[1].tolist():
                    lst.append(v2[0])
            s.add_clause(lst)
            self.num_clauses += 1

        self.times[6] = time()

        # if a transform is used in an interior corona, its halo cells must
        # be used
        for k, v in self.transforms.items():

            # do not consider nth corona
            if k[0] == self.k_cor:
                continue

            # consider all cells in halo of transform
            for c in v[2]:
                s.add_clause([-v[0], self.cells[tuple(c)]])
                self.num_clauses += 1

        self.times[7] = time()

        # print(cpu_count())
        # print('here0')
        # with Manager() as manager:
        #     lst = manager.list()
        with Pool(processes=cpu_count()) as pool:
            for lst in pool.imap_unordered(self.check_overlap_mp, range(0, self.grid.size[0])):
                for l in lst:
                    s.add_clause(l)
                    self.num_clauses += 1
        # print('here1.5')

        # print(time() - ts)

        # TODO - optimize (takes a long time ~2000s for k=2)
        #  probably caused by the is_overlapping method
        # used transforms cannot overlap
        # for k1, v1 in self.transforms.items():
        #     for k2, v2 in self.transforms.items():
        #
        #         # only consider pairings where idx1 < idx2 so unnecessary
        #         # clauses are not generated
        #         if v1[0] >= v2[0]:
        #             continue
        #
        #         # checks if any two rows in the transform are the same
        #         if self.grid.is_overlapping(v1[1], v2[1]):
        #             s.add_clause([-v1[0], -v2[0]])

        # same speed as above
        # for (k1, v1), (k2, v2) in itertools.combinations(self.transforms.items(), 2):
        #     # if np.max(np.abs(v1[1][0] - v2[1][0])) > 2*self.shape_rad:
        #     # continue
        #
        #     # seems to not change speed at k=0,1?
        #     if max(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * self.shape_rad:
        #         continue
        #
        #     # checks if any two rows in the transform are the same
        #     if self.grid.is_overlapping(v1[1], v2[1]):
        #         s.add_clause([-v1[0], -v2[0]])

        self.times[8] = time()

        # TODO - optimize (takes a long time ~30s for k=2)
        #  probably caused by the is_overlapping method
        # if a transform is used in a k corona, it must be adjacent to one
        # in a k-1 corona
        # for k1, v1 in self.transforms.items():
        #     if k1[0] == 0:
        #         continue
        #
        #     lst = [-v1[0]]
        #     for k2, v2 in self.transforms.items():
        #
        #         # only consider the k and k-1 coronas
        #         if k1[0] != k2[0] + 1:
        #             continue
        #
        #         # if np.max(np.abs(v1[1][0] - v2[1][0])) > 2 * (shape_rad + 1):
        #         #     continue
        #
        #         # checks if any two rows in the transform and second transform's
        #         # halo are the same
        #         if self.grid.is_overlapping(v1[1], v2[2]):
        #             lst.append(v2[0])
        # print('here2')
        # with Manager() as manager:
        #     lst = manager.list()
        with Pool(processes=cpu_count()) as pool:
            for lst in pool.imap_unordered(self.check_halo_overlap_mp, range(0, self.grid.size[0])):
                for l in lst:
                    # print(l)
                    s.add_clause(l)
                    self.num_clauses += 1
        # print('here3')

        self.times[9] = time()

        # if a transform is used in a k corona, it cannot be adjacent to one
        # in a 0...k-2 corona
        for k1, v1 in self.transforms.items():
            for k2, v2 in self.transforms.items():

                # only consider coronas < k-1
                if k2[0] >= k1[0] - 1:
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):
                    s.add_clause([-v1[0], -v2[0]])
                    self.num_clauses += 1

        self.times[10] = time()
        #
        # # TODO - add hole suppression
        # # not sure why we can't just have a and b and c and d => e
        # # where a, b, c, d, are the tiles adjacent to e
        for k, v in self.cells.items():
            lst = []

            # find all adjacent cells
            for c in self.grid.haloIdx(k):
                lst.append(-self.cells[tuple(c)])

            lst.append(v)
            s.add_clause(lst)
            self.num_clauses += 1

        self.times[11] = time()

        self.sat = s
        return self.sat

    def solve_sat(self) -> list:
        self.times[12] = time()

        if self.sat.solve():
            self.model = self.sat.get_model()
        else:
            self.model = None

        self.times[13] = time()

        return self.model

    def clear_sat(self):
        self.sat.delete()
        self.sat = None
        self.model = None

    def get_transform(self, idx: int):

        idx -= 1
        if idx >= (self.k_cor + 1) * self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices):
            return None

        v1 = int(idx / (self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)))
        idx -= v1 * (self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices))
        v2 = int(idx / (self.grid.size[1] * len(self.rotation_matrices)))
        idx -= v2 * (self.grid.size[1] * len(self.rotation_matrices))
        v3 = int(idx / (len(self.rotation_matrices)))
        idx -= v3 * (len(self.rotation_matrices))

        return v1, v2, v3, idx

    def get_transform_idx(self, key: tuple[int, int, int, int]):
        return key[0] * self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices) + \
               key[1] * self.grid.size[1] * len(self.rotation_matrices) + \
               key[2] * len(self.rotation_matrices) + \
               key[3] + \
               1

    def check_overlap_mp(self, i1: int):
        out = []
        for (k1, v1), (k2, v2) in itertools.combinations(self.transforms.items(), 2):
            if k1[1] != i1 and k2[1] != i1:
                continue

            # seems to not change speed at k=0,1?
            if max(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * self.shape_rad:
                continue

            # checks if any two rows in the transform are the same
            if self.grid.is_overlapping(v1[1], v2[1]):
                out.append([-v1[0], -v2[0]])

        return out

    def check_halo_overlap_mp(self, i1: int):
        out = []
        for k1, v1 in self.transforms.items():
            if k1[1] != i1:
                continue
            if k1[0] == 0:
                continue

            l = [-v1[0]]
            for k2, v2 in self.transforms.items():

                # only consider the k and k-1 coronas
                if k1[0] != k2[0] + 1:
                    continue

                if max(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * (self.shape_rad + 1):
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):
                    l.append(v2[0])

            if len(l) > 0:
                out.append(l)

        return out

    def check_rotational_symmetry(self):

        seen = []
        out = [1] * len(self.rotation_matrices)

        for idx, rm in enumerate(self.rotation_matrices):
            t = set(map(tuple, np.matmul(rm, self.shape.T).T))
            if t in seen:
                out[idx] = 0
            else:
                seen.append(t)

        return out

    def write(self, directory='.', plot=True):
        filename = str(int(time()))

        with open(directory + '/' + filename + '.txt', 'w+') as f:

            f.write(f'File: {filename}\n')
            f.write(f'Start Time: {self.times[0]}\n')
            f.write(f'Grid Size: {self.grid.size}\n')
            f.write(f'Coronas: {self.k_cor}\n')
            f.write(f'Shape Size: {self.shape_size}\n')
            f.write(f'Shape Radius: {self.shape_rad}\n')
            f.write(f'Shape: \n')
            f.write(f'{str(self.shape)}\n')

            f.write('--------------------------------------------------\n')

            f.write(f'Total Clauses: {self.num_clauses}\n')
            f.write(f'SAT Model: \n')
            if self.model is not None:
                f.write(f'\t{str(self.model)}\n')
            else:
                f.write(f'\tNo Model (Unsolvable)\n')
            f.write('Times: \n')
            sum_time = 0.0
            for i, t in enumerate(self.times):
                f.write(f'\tTimestamp {i:02d}: {t - sum_time - self.times[0]}\n')
                sum_time = t - self.times[0]
            f.write(f'Total Time: {sum_time}\n')

            f.write('--------------------------------------------------\n')

            f.write('Transforms: \n')
            if self.model is not None:
                for i in self.model:
                    if i <= 0:
                        continue
                    t = self.get_transform(i)
                    if t is None:
                        break
                    f.write(f"Transform ID: {i} {self.get_transform(i)}: \n")
                    f.write(f'{str(self.transforms[t][1])}\n')
            else:
                f.write('\tNo Model (Unsolvable)')

        if plot:
            self.plot(show=False, write=True, filename=filename, directory=directory)

    @abstractmethod
    def plot(self, show, write, filename, directory):
        pass
