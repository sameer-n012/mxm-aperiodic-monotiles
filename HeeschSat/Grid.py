import numpy as np
from abc import ABC, abstractmethod
from typing import Union, Optional


class Grid(ABC):

    def __init__(self, size: tuple, basis: np.array):
        """
        Creates a grid structure based on the size provided. The size provided
        should be in the basis given. The basis should be a 2-dimensional array
        as such: [
            [ u_1, u_1 ],
            [ v_2, v_2 ]
        ] where each row is a basis vector

        """
        self.basis = np.array([1, 0], [0, 1]) if basis is None else basis.T
        self.size = size
        self.grid = np.full(self.size, False)  # TODO - delete unused

    @staticmethod
    @abstractmethod
    def is_adjacent(x: Union[tuple, list, np.array], y: Union[tuple, list, np.array]):
        pass

    # TODO - delete unused
    def halo(self, cell: tuple[int]) -> np.array:
        """
        Returns the cells that are adjacent to the current cell based on the
        adjacency function the grid was created with. The cells returned will
        not be in any particular order.
        """

        return np.array([self.grid[tuple(x)] for x in np.transpose(np.nonzero(self.haloIdx(cell)))])

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
        ] where each row is a distinct cell
        """

        arr = np.array([self.is_adjacent(i, cell) for i in self.indices()]).reshape(self.size)
        return np.transpose(np.nonzero(arr))

    def apply_basis(self, coords: np.ndarray) -> np.array:
        """
        Takes the coordinates in the grid's specified basis and converts it to
        coordinates in the standard basis for visualization. The coordinates
        should be in a 2-dimensional array as such: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ] where each row is a distinct point

        """
        return np.matmul(self.basis, coords.T)

    def indices(self) -> np.array:
        """
        Returns an array of indices on the grid
        """

        return np.array([(i, j) for j in range(self.size[1]) for i in range(self.size[0])])

    def is_overlapping(self, m1: np.ndarray, m2: np.ndarray) -> bool:
        """
        Checks whether the two matrices, m1 and m2 have any similar rows. Used
        to test if there are overlapping points in two shapes. It is required
        that the matrices have no duplicate rows to begin withThe matrices
        should be in the form: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ]
        """

        stack = np.concatenate((m1, m2), axis=0)
        _, c = np.unique(stack, axis=0, return_counts=True)
        return (c != 1).all()

    def in_bounds(self, shape: np.ndarray) -> bool:
        """
        Checks whether the shape given is within the bounds of the grid. The
        shape given should be of the form: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ] where each row is a cell the shape occupies
        """

        return ((shape[:, 0] < self.size[0]).all() and
                (shape[:, 1] < self.size[0]).all() and
                (shape > 0).all())
