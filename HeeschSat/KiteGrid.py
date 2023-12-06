import numpy as np

from HeeschSat.Grid import Grid
from math import sqrt

from collections import deque


class KiteGrid(Grid):

    def __init__(self, size: tuple, pos: int):
        super().__init__(size, None)
        self.pos = pos # TODO unused?

    def apply_basis(self, coords: np.ndarray) -> np.array:
        """
        Turn n sets of 3-dim-coords into coords in standard basis R^2,
        will be given as an n*3 matrix
        """
        # TODO complete
        pass

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
            new_kite.append((kite[2] + 1) % 6)
            transformation.append(new_kite)

        return transformation

    @staticmethod
    def _reflect(tile: list):
        ref_mat = np.array([[-1, 0],
                            [0, 1]])
        # transformation will be a list of each original kite reflected.
        transformation = []
        for kite in tile:
            # multiply first two coords by reflection matrix
            new_kite = np.dot(ref_mat, kite[:2])
            # handles 3rd coordinate
            match kite[2]:
                case 0:
                    new_kite.append(2)
                case 1:
                    new_kite.append(1)
                case 2:
                    new_kite.append(0)
                case 3:
                    new_kite.append(5)
                case 4:
                    new_kite.append(4)
                case 5:
                    new_kite.append(3)
            transformation.append(new_kite)

        return transformation

    @staticmethod
    def get_transformations(tile: list):
        # will contain each possible transformation of given tile
        stack = deque()
        # add original tile and reflected tile to stack
        stack.append(KiteGrid._reflect(tile))
        stack.append(tile)
        while KiteGrid._rotate(stack[-1]) != tile:
            # find reflected & non-reflected rotation of most recently added tile
            ref_rot_tile = KiteGrid._reflect(KiteGrid._rotate(stack[-1]))
            rot_tile = KiteGrid._rotate(stack[-1])
            # push new tiles to stack
            stack.append(ref_rot_tile)
            stack.append(rot_tile)

        return list(stack)
