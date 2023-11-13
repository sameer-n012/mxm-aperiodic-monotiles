import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt


class KiteGrid(Grid):

    def __init__(self, size: tuple, pos: int):
        super().__init__(size, np.array([[1, 0], [0.5, sqrt(3) * 0.5]]))
        self.pos = pos

    @staticmethod
    def is_adjacent(x, y):
        if x[0] == y[0] and x[1] == y[1]:
            return False
        match x[2]:
            case 0:
                return (
                        x[0] + 1 == y[0] and
                        ((x[1] == y[1] and (y[2] == 3 or y[2] == 4)) or
                         (x[1] == y[1] - 1 and (y[2] == 2 or y[2] == 3)))
                )
            case 1:
                return (
                        (x[0] == y[0] and x[1] + 1 == y[1] and (y[2] == 5 or y[2] == 4)) or
                        (x[0] + 1 == y[0] and x[1] == y[1] and (y[2] == 3 or y[2] == 4))
                )
            case 2:
                return (
                        x[1] + 1 == y[1] and
                        ((x[0] == y[0] and (y[2] == 4 or y[2] == 5)) or
                         (x[0] - 1 == y[0] and (y[2] == 5 or y[2] == 0)))
                )
            case 3:
                return (
                    x[0] - 1 == y[0] and
                    ((x[1] == y[1] and (y[2] == 1 or y[2] == 0)) or
                     (x[1] + 1 == y[1] and (y[2] == 5 or y[2] == 0)))
                )
            case 4:
                return (
                    (x[0] - 1 == y[0] and x[1] == y[1] and (y[2] == 1 or y[2] == 0)) or
                    (x[0] == y[0] and x[1] - 1 == y[1] and (y[2] == 1 or y[2] == 2))
                )
            case 5:
                return (
                    x[1] - 1 == y[1] and
                    ((x[0] == y[0] and (y[2] == 1 or y[2] == 2)) or
                     (x[0] + 1 == y[0] and (y[2] == 2 or y[2] == 3)))
                )
