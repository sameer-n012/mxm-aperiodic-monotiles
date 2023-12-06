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
        ] where each row is a basis vector. Runs in O(1) time.

        """

        self.basis = np.array([1, 0], [0, 1]) if basis is None else basis.T
        self.size = size

    @staticmethod
    @abstractmethod
    def is_adjacent(x: Union[tuple, list, np.array], y: Union[tuple, list, np.array]):
        pass

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

        arr = np.array([self.is_adjacent(i, cell) for i in self.indices()]).reshape(self.size, order='F')
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
        ] where each row is a distinct point. Runs in O(mn^2) time for a grid
        of size n*n and a coordinate shape of size n*2.
        """

        return np.matmul(self.basis, coords.T)

    def indices(self) -> np.array:
        """
        Returns an array of indices on the grid. Runs in O(n^2) time for a grid
        of size n*n.
        """

        return np.array([(i, j) for j in range(self.size[1]) for i in range(self.size[0])])

    def is_overlapping(self, m1: np.ndarray, m2: np.ndarray) -> bool:
        """
        Checks whether the two matrices, m1 and m2 have any similar rows. Used
        to test if there are overlapping points in two shapes. It is required
        that the matrices have no duplicate rows to begin with. The matrices
        should be in the form: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ]. Runs in O(m) time for matrices of size m*2.
        """

        # stack = np.concatenate((m1, m2), axis=0)
        # _, c = np.unique(stack, axis=0, return_counts=True)
        # return not (c == 1).all()
        return len(set(map(tuple, m1)) & set(map(tuple, m2))) > 0

    def is_overlapping_set(self, m1: set, m2: set) -> bool:
        return len(m1 & m2) > 0

    def in_bounds(self, shape: np.ndarray) -> bool:
        """
        Checks whether the shape given is within the bounds of the grid. The
        shape given should be of the form: [
            [ x_1, y_1 ],
            [ x_2, y_2 ],
            ...
            [ x_n, y_n ]
        ] where each row is a cell the shape occupies. Runs in O(m) time for
        a shape of size m*2.
        """

        return ((shape[:, 0] < self.size[0]).all() and
                (shape[:, 1] < self.size[1]).all() and
                (shape >= 0).all())
