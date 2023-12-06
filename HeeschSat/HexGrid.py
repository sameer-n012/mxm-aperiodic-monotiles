import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt


class HexGrid(Grid):

    def __init__(self, size: tuple):
        super().__init__(size, np.array([[1, 0], [0.5, sqrt(3) * 0.5]]))

    def indices(self) -> np.ndarray:
        """
        Returns an array of indices on the grid. Runs in O(n^2) time for a grid
        of size n*n.
        """
        return np.array([(i, j) for j in range(self.size[1]) for i in range(self.size[0])])

    @staticmethod
    def is_adjacent(x, y) -> bool:
        """
        Determines whether two coordinates are adjacent to each other
        on the grid. The coordinates x, y should be specified as a tuple
        or list of length 2: [x0, x1]. Two coordinates are adjacent if any of
        the following are true:
            - x0 == y0 and x1 and y1 differ by at most 1
            - x1 == y1 and x0 and y0 differ by at most 1
            - x0 is 1 greater than y0 and x1 is 1 less than y1
            - x0 is 1 less than y0 and x1 is 1 greater than y1
        Runs in O(1) time. Returns true if they are adjacent and false
        otherwise.
        """
        return (
                (x[0] == y[0] and x[1] - 1 <= y[1] <= x[1] + 1) or
                (x[1] == y[1] and x[0] - 1 <= y[0] <= x[0] + 1) or
                (x[0] - 1 == y[0] and x[1] + 1 == y[1]) or
                (x[0] + 1 == y[0] and x[1] - 1 == y[1])
        )

    @staticmethod
    def get_adjacent(x) -> np.ndarray:
        pass


