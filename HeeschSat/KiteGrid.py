import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt

from collections import deque


class KiteGrid(Grid):

    def __init__(self, size: tuple, pos: int):
        super().__init__(size, np.array([[1, 0, 0], [0.5, sqrt(3) * 0.5, 0], [0, 0, 1]]))
        self.pos = pos  # TODO unused?

    # def apply_basis(self, coords: np.ndarray) -> np.array:
    #     """
    #     Turn n sets of 3-dim-coords into coords in standard basis R^2,
    #     will be given as an n*3 matrix
    #     """
    #     # TODO complete
    #     pass

    def haloIdx(self, cell: tuple[int]) -> np.ndarray:
        """
        Returns the indices of the cells that are adjacent to the current cell
        based on the adjacency function the grid was created with. The cells
        returned will not be in any particular order. The cells returned will
        be in the form: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ] where each row is a distinct cell. Runs in O(n^2) time for a grid of
        size n*n.
        """

        out = []
        for j in range(self.size[1]):
            for i in range(self.size[0]):
                for k in range(self.size[2]):
                    if self.is_adjacent([i, j, k], cell):
                        out.append([i, j, k])

        # arr = np.array([self.is_adjacent(i, cell) for i in self.indices()]).reshape(self.size, order='C')
        # return np.transpose(np.nonzero(arr))
        return out

    def indices(self) -> np.ndarray:
        """
        Returns an array of indices on the grid. Runs in O(kn^2) time for a grid
        of size n*n*k.
        """
        return np.array([(i, j, k) for j in range(self.size[1]) for i in range(self.size[0])
                         for k in range(self.size[2])])

    @staticmethod
    def is_adjacent(x, y):

        # TODO - do we count corners touching as adjacent??

        # case where two kites are in same hex
        if x[0] == y[0] and x[1] == y[1]:
            # return np.abs(x[2] - y[2]) == 1 or np.abs(x[2] - y[2]) == 5
            return True

        match x[2]:
            case 0:
                return (
                    (x[0] + 1 == y[0] and x[1] == y[1] and (y[2] == 2 or y[2] == 3)) or
                    (x[0] == y[0] and x[1] + 1 == y[1] and (y[2] == 3 or y[2] == 4))
                )
            case 1:
                return (
                    (x[0] == y[0] and x[1] + 1 == y[1] and (y[2] == 3 or y[2] == 4)) or
                    (x[0] - 1 == y[0] and x[1] + 1 == y[1] and (y[2] == 4 or y[2] == 5))
                )
            case 2:
                return (
                    (x[0] - 1 == y[0] and x[1] + 1 == y[1] and (y[2] == 4 or y[2] == 5)) or
                    (x[0] - 1 == y[0] and x[1] == y[1] and (y[2] == 0 or y[2] == 5))
                )
            case 3:
                return (
                    (x[0] - 1 == y[0] and x[1] == y[1] and (y[2] == 0 or y[2] == 5)) or
                    (x[0] == y[0] and x[1] - 1 == y[1] and (y[2] == 0 or y[2] == 1))
                )
            case 4:
                return (
                    (x[0] == y[0] and x[1] - 1 == y[1] and (y[2] == 0 or y[2] == 1)) or
                    (x[0] + 1 == y[0] and x[1] - 1 == y[1] and (y[2] == 1 or y[2] == 2))
                )
            case 5:
                return (
                    (x[0] + 1 == y[0] and x[1] - 1 == y[1] and (y[2] == 1 or y[2] == 2)) or
                    (x[0] + 1 == y[0] and x[1] == y[1] and (y[2] == 2 or y[2] == 3))
                )
            case _:
                return False

        # match x[2]:
        #     case 0:
        #         return (
        #                 x[0] + 1 == y[0] and
        #                 ((x[1] == y[1] and (y[2] == 3 or y[2] == 4)) or
        #                  (x[1] == y[1] - 1 and (y[2] == 2 or y[2] == 3)))
        #         )
        #     case 1:
        #         return (
        #                 (x[0] == y[0] and x[1] + 1 == y[1] and (y[2] == 5 or y[2] == 4)) or
        #                 (x[0] + 1 == y[0] and x[1] == y[1] and (y[2] == 3 or y[2] == 4))
        #         )
        #     case 2:
        #         return (
        #                 x[1] + 1 == y[1] and
        #                 ((x[0] == y[0] and (y[2] == 4 or y[2] == 5)) or
        #                  (x[0] - 1 == y[0] and (y[2] == 5 or y[2] == 0)))
        #         )
        #     case 3:
        #         return (
        #                 x[0] - 1 == y[0] and
        #                 ((x[1] == y[1] and (y[2] == 1 or y[2] == 0)) or
        #                  (x[1] + 1 == y[1] and (y[2] == 5 or y[2] == 0)))
        #         )
        #     case 4:
        #         return (
        #                 (x[0] - 1 == y[0] and x[1] == y[1] and (y[2] == 1 or y[2] == 0)) or
        #                 (x[0] == y[0] and x[1] - 1 == y[1] and (y[2] == 1 or y[2] == 2))
        #         )
        #     case 5:
        #         return (
        #                 x[1] - 1 == y[1] and
        #                 ((x[0] == y[0] and (y[2] == 1 or y[2] == 2)) or
        #                  (x[0] + 1 == y[0] and (y[2] == 2 or y[2] == 3)))
        #         )

    @staticmethod
    def _rotate(tile: list):
        rot_mat = np.array([[0, -1],
                            [1, 1]])
        # transformation will be an array of each original kite rotated once.
        transformation = []
        for kite in tile:
            # multiply first two coords by rotation matrix
            new_kite = np.dot(rot_mat, kite[:2])
            # handles 3rd coordinate
            # new_kite.append((kite[2] + 1) % 6)
            new_kite = np.hstack((new_kite, [(kite[2] + 1) % 6]))
            transformation.append(new_kite)

        return transformation

    @staticmethod
    def _reflect(tile: list):
        ref_mat = np.array([[-1, 0],
                            [1, 1]])
        # transformation will be a list of each original kite reflected.
        transformation = []
        for kite in tile:
            # multiply first two coords by reflection matrix
            new_kite = np.dot(ref_mat, kite[:2])
            # handles 3rd coordinate
            match kite[2]:
                case 0:
                    # new_kite.append(2)
                    new_kite = np.hstack((new_kite, [1]))
                case 1:
                    # new_kite.append(1)
                    new_kite = np.hstack((new_kite, [0]))
                case 2:
                    # new_kite.append(0)
                    new_kite = np.hstack((new_kite, [5]))
                case 3:
                    # new_kite.append(5)
                    new_kite = np.hstack((new_kite, [4]))
                case 4:
                    # new_kite.append(4)
                    new_kite = np.hstack((new_kite, [3]))
                case 5:
                    # new_kite.append(3)
                    new_kite = np.hstack((new_kite, [2]))
            transformation.append(new_kite)

        return transformation

    @staticmethod
    def get_transformations(tile: list):
        # will contain each possible transformation of given tile
        stack = deque()
        # add original tile and reflected tile to stack
        stack.append(KiteGrid._reflect(tile))
        stack.append(tile)
        while (KiteGrid._rotate(stack[-1]) != tile).any():
            # find reflected & non-reflected rotation of most recently added tile
            ref_rot_tile = KiteGrid._reflect(KiteGrid._rotate(stack[-1]))
            rot_tile = KiteGrid._rotate(stack[-1])
            # push new tiles to stack
            stack.append(ref_rot_tile)
            stack.append(rot_tile)

        return list(stack)
