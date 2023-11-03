import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt


class KiteGrid(Grid):

    def __init__(self, size: tuple, pos: int):
        super().__init__(size, np.array([[1, 0], [0.5, sqrt(3) * 0.5]]))
        self.pos = pos

    @staticmethod
    # Called by is_adjacent if pos of x is 0
    def is_adj_pos0(x, y):
        return (
                x[0] + 1 == y[0] and
                ((x[1] == y[1] and (y[2] == 3 or y[2] == 4)) or
                 (x[1] == y[1] - 1 and (y[2] == 2 or y[2] == 3)))
        )

# TODO: Write helper methods 1 through 5
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

    @staticmethod
    def is_adjacent(x, y):
        match x[2]:
            case 0:
                return KiteGrid.is_adj_pos0(x, y)
            case 1:
                return KiteGrid.is_adj_pos1(x, y)
            case 2:
                return KiteGrid.is_adj_pos2(x, y)
            case 3:
                return KiteGrid.is_adj_pos3(x, y)
            case 4:
                return KiteGrid.is_adj_pos4(x, y)
            case 5:
                return KiteGrid.is_adj_pos5(x, y)
            case _:
                raise Exception("Position of kite may not be 6.")
