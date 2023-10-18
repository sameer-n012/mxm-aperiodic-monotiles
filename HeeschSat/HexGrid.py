import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt


class HexGrid(Grid):

    def __init__(self, size: tuple):
        super().__init__(size, np.array([[1, 0], [0.5, sqrt(3) * 0.5]]))

    @staticmethod
    def is_adjacent(x, y):
        return (
                (x[0] == y[0] and x[1] - 1 <= y[1] <= x[1] + 1) or
                (x[1] == y[1] and x[0] - 1 <= y[0] <= x[0] + 1) or
                (x[0] - 1 == y[0] and x[1] + 1 == y[1]) or
                (x[0] + 1 == y[0] and x[1] - 1 == y[1])
        )


