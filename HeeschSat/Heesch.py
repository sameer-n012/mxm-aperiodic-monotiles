from abc import ABC, abstractmethod
import numpy as np
import itertools
from pysat.solvers import Solver

class Heesch(ABC):

    def __init__(self):
        self.k_cor = None
        self.grid = None
        self.transforms = {}
        self.cells = {}
        self.rotation_matrices = []
        self.shape = None
        self.shape_size = None

    def generate_variables(self):

        # transform variables (if a transform is used)
        for idx, val in enumerate(itertools.product(range(0, self.k_cor),
                                                    range(0, self.grid.size[0]),
                                                    range(0, self.grid.size[1]),
                                                    range(0, len(self.rotation_matrices))
                                                    )):
            translate_mat = np.column_stack((np.full(np.shape_size, val[1]), np.full(np.shape_size, val[2])))
            rotate_mat = self.rotation_matrices[val[3]]
            transform = np.matmul(rotate_mat, self.shape.T).T + translate_mat
            halo = np.unique([i for c in self.shape for i in self.grid.haloIdx(c)])
            self.transforms[val] = idx, transform, halo

        # cell variables (if a cell is taken)
        offset = len(self.transforms)
        for idx, i in enumerate(self.grid.indices()):
            self.cells[i] = offset + idx

    @abstractmethod
    def construct_sat(self) -> Solver:
        """
        Note that a -> b | c === -a | b | c
        """

        s = Solver(name='g3')

        # 0-corona always used
        k_0 = self.transforms[(0, self.grid.size[0]/2, self.grid.size[1]/2, 0)]
        s.add_clause([k_0])

        # if a transform is used, its cells are used
        for k, v in self.transforms:
            for c in v[2]:
                s.add_clause([-v[0], self.cells[c]])

        # if a cell is used, some transform must use it


        c1 = [-1, 2]
        s.add_clause(c1)

        return s

    @abstractmethod
    def solve_sat(self, sat: Solver) -> list:
        if sat.solve():
            m = sat.get_model()
        else:
            m = None

        sat.delete()
        return m










