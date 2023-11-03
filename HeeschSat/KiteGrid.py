import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt


class KiteGrid(Grid):

    def __init__(self, size: tuple, pos: int):
        super().__init__(size, np.array([[1, 0], [0.5, sqrt(3) * 0.5]]))
        self.pos = pos

    # TODO: Clean up formatting
    @staticmethod
    def is_adj_pos1(x, y):
        return (
            None
        )

    @staticmethod
    def is_adj_pos2(x, y):
        return (
            None
        )

    @staticmethod
    def is_adj_pos3(x, y):
        return (
            None
        )

    @staticmethod
    def is_adj_pos4(x, y):
        return (
            None
        )

    @staticmethod
    def is_adj_pos5(x, y):
        return (
            None
        )

# TODO: finish cases 3, 4, 5
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
                        (x[0] == y[0] and (y[2] == 4 or y[2] == 5)) or
                        (x[0] - 1 == y[0] and (y[2] == 5 or y[2] == 0))
                )
            case 3:
                return None
            case 4:
                return None
            case 5:
                return None
            case _:
                return False
