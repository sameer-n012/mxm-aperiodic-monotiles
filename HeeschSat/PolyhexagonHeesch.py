from HeeschSat.Heesch import Heesch
from HeeschSat.HexGrid import HexGrid
import numpy as np
from numpy import identity
import itertools


class PolyhexagonHeesch(Heesch):

    def __init__(self, shape, coronas, grid_size):
        super().__init__(coronas)
        self.grid = HexGrid(grid_size)
        self.rotation_matrices = [
            np.array([[1, 0],
                      [0, 1]])
        ]
        self.k_cor = coronas
        self.shape = shape
        self.shape_size = self.shape.shape
