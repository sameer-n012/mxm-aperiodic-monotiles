from abc import ABC, abstractmethod
import numpy as np
import itertools
from pysat.solvers import Solver


class Heesch(ABC):

    def __init__(self, coronas):
        self.k_cor = coronas
        self.grid = None
        self.transforms = {}
        self.cells = {}
        self.rotation_matrices = []
        self.shape = None
        self.shape_size = None
        self.sat = Solver(name="g3")
        self.model = None

    def generate_variables(self):

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
                self.transforms[tuple(val)] = idx, transform, halo

        # cell variables (if a cell is taken)
        offset = len(self.transforms)
        for idx, i in enumerate(self.grid.indices()):
            self.cells[tuple(i)] = offset + idx

    def construct_sat(self) -> Solver:
        """
        Note that a, b => c | d === -a | -b | c | d
        """

        s = self.sat

        # 0-corona always used
        k_0 = self.transforms[(0, self.grid.size[0] / 2, self.grid.size[1] / 2, 0)]
        s.add_clause([k_0])

        # if a transform is used, its cells are used
        for k, v in self.transforms:

            # ad each cell in the transform
            for c in v[1]:
                s.add_clause([-v[0], self.cells[c]])

        # if a cell is used, some transform must use it
        for k1, v1 in self.cells:
            lst = [-v1[0]]

            # check if the cell is in the transform for
            # each transform
            for k2, v2 in self.transforms:
                if k1 in v2[1]:
                    lst.append(v2[0])
            s.add_clause(lst)

        # if a transform is used in an interior corona, its halo cells must
        # be used
        for k, v in self.transforms:

            # do not consider nth corona
            if k[0] == self.k_cor:
                continue

            # consider all cells in halo of transform
            for c in v[2]:
                s.add_clause([-v[0], self.cells[c]])

        # used transforms cannot overlap
        for k1, v1 in self.transforms:
            for k2, v2 in self.transforms:

                # only consider pairings where idx1 < idx2 so unnecessary
                # clauses are not generated
                if v1[0] >= v2[0]:
                    continue

                # checks if any two rows in the transform are the same
                if self.grid.is_overlapping(v1[1], v2[1]):
                    s.add_clause([-v1[0], -v2[0]])

        # if a transform is used in a k corona, it must be adjacent to one
        # in a k-1 corona
        for k1, v1 in self.transforms:
            lst = [-v1[0]]
            for k2, v2 in self.transforms:

                # only consider the k and k-1 coronas
                if k1[0] != k2[0] + 1:
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):
                    lst.append(v2[0])
            s.add_clause(lst)

        # if a transform is used in a k corona, it cannot be adjacent to one
        # in a 0...k-2 corona
        for k1, v1 in self.transforms:
            for k2, v2 in self.transforms:

                # only consider coronas < k-1
                if k2[0] >= k1[0] - 1:
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):
                    s.add_clause([-v1[0], -v2[0]])

        # TODO - add hole suppression
        # not sure why we can't just have a and b and c and d => e
        # where a, b, c, d, are the tiles adjacent to e
        for k, v in self.cells:
            lst = []

            # find all adjacent cells
            for c in self.grid.haloIdx(k):
                lst.append(-self.cells[c])

            lst.append(v[0])
            s.add_clause(lst)

        return s

    def solve_sat(self) -> list:
        if self.sat.solve():
            self.model = self.sat.get_model()
        else:
            self.model = None

        return self.model

    def clear_sat(self):
        self.sat.delete()
        self.sat = None
        self.model = None
