from abc import ABC, abstractmethod
import numpy as np
import itertools
from pysat.solvers import Solver
from time import time


class Heesch(ABC):

    def __init__(self, coronas):
        self.k_cor = coronas
        self.grid = None
        self.transforms = {}
        self.cells = {}
        self.rotation_matrices = []
        self.shape = None
        self.shape_size = None
        self.sat = None
        self.model = None

    def generate_variables(self):

        max_transforms = (self.k_cor + 1)*self.grid.size[0]*self.grid.size[1]*len(self.rotation_matrices)

        # transform variables (if a transform is used)
        for idx, val in enumerate(itertools.product(range(0, self.k_cor + 1),
                                                    range(0, self.grid.size[0]),
                                                    range(0, self.grid.size[1]),
                                                    range(0, len(self.rotation_matrices))
                                                    )):
            translate_mat = np.empty(self.shape_size, dtype=int)
            translate_mat[:, 0] = val[1]
            translate_mat[:, 1] = val[2]
            rotate_mat = self.rotation_matrices[val[3]]
            transform = np.matmul(rotate_mat, self.shape.T).T + translate_mat

            # only include transforms that fall within the bounds of the grid
            if self.grid.in_bounds(transform):
                halo = np.unique([i for c in self.shape for i in self.grid.haloIdx(c)], axis=0)
                self.transforms[tuple(val)] = idx + 1, transform, halo

        # cell variables (if a cell is taken)
        offset = max_transforms + 1
        for idx, i in enumerate(self.grid.indices()):
            self.cells[tuple(i)] = offset + idx

    def construct_sat(self) -> Solver:
        """
        Note that a, b => c | d === -a | -b | c | d
        """

        s = Solver(name="g3")

        ts = time()

        # 0-corona always used
        k_0 = self.transforms[(0, int(self.grid.size[0] / 2), int(self.grid.size[1] / 2), 0)]
        s.add_clause([k_0[0]])

        print(time() - ts)

        # if a transform is used, its cells are used
        for k, v in self.transforms.items():

            # ad each cell in the transform
            for c in v[1]:
                s.add_clause([-v[0], self.cells[tuple(c)]])

        print(time() - ts)

        # if a cell is used, some transform must use it
        for k1, v1 in self.cells.items():
            lst = [-v1]

            # check if the cell is in the transform for
            # each transform
            for k2, v2 in self.transforms.items():
                if list(k1) in v2[1].tolist():
                    lst.append(v2[0])
            s.add_clause(lst)

        print(time() - ts)


        # if a transform is used in an interior corona, its halo cells must
        # be used
        for k, v in self.transforms.items():

            # do not consider nth corona
            if k[0] == self.k_cor:
                continue

            # consider all cells in halo of transform
            for c in v[2]:
                s.add_clause([-v[0], self.cells[tuple(c)]])

        print(time() - ts)

        # TODO - optimize (takes a long time)
        # used transforms cannot overlap
        for k1, v1 in self.transforms.items():
            for k2, v2 in self.transforms.items():

                # only consider pairings where idx1 < idx2 so unnecessary
                # clauses are not generated
                if v1[0] >= v2[0]:
                    continue

                # checks if any two rows in the transform are the same
                if self.grid.is_overlapping(v1[1], v2[1]):
                    s.add_clause([-v1[0], -v2[0]])

        print(time() - ts)

        # if a transform is used in a k corona, it must be adjacent to one
        # in a k-1 corona
        for k1, v1 in self.transforms.items():
            lst = [-v1[0]]
            for k2, v2 in self.transforms.items():

                # only consider the k and k-1 coronas
                if k1[0] != k2[0] + 1:
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):
                    lst.append(v2[0])
            if len(lst) > 1:
                s.add_clause(lst)

        print(time() - ts)

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

        print(time() - ts)
        #
        # # TODO - add hole suppression
        # # not sure why we can't just have a and b and c and d => e
        # # where a, b, c, d, are the tiles adjacent to e
        # for k, v in self.cells.items():
        #     lst = []
        #
        #     # find all adjacent cells
        #     for c in self.grid.haloIdx(k):
        #         lst.append(-self.cells[tuple(c)])
        #
        #     lst.append(v)
        #     s.add_clause(lst)

        print(time() - ts)

        self.sat = s
        return self.sat

    def solve_sat(self) -> list:
        if self.sat.solve():
            print('solved sat')
            self.model = self.sat.get_model()
        else:
            print('unsolvable')
            self.model = None

        return self.model

    def clear_sat(self):
        self.sat.delete()
        self.sat = None
        self.model = None

    def get_transform(self, idx):

        idx -= 1
        v1 = int(idx / (self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)))
        idx -= v1*(self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices))
        v2 = int(idx / (self.grid.size[1] * len(self.rotation_matrices)))
        idx -= v2*(self.grid.size[1] * len(self.rotation_matrices))
        v3 = int(idx / (len(self.rotation_matrices)))
        idx -= v3*(len(self.rotation_matrices))

        return v1, v2, v3, idx
