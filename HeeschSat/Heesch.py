from abc import ABC, abstractmethod
import numpy as np
import itertools
from pysat.solvers import Solver
from time import time
import matplotlib.colors as mcolors
from multiprocessing import Pool, cpu_count

"""
Class implemented based on Kaplan's Heesch Numbers of Unmarked Polyforms
(https://cs.uwaterloo.ca/~csk/heesch/unmarked.pdf)
"""


# defines the colors for plotting
def define_colors():
    lst = sorted(
        mcolors.CSS4_COLORS, key=lambda c: tuple(mcolors.rgb_to_hsv(mcolors.to_rgb(c)))
    )
    lst.remove("black")
    return lst, [0, 14, 33, 53, 64, 85, 122, 132]


class Heesch(ABC):
    # plot_colors = ['k', 'b', 'r', 'y', 'g', 'c', 'm']
    plot_colors, c_majors = define_colors()
    c_spacing = 3
    num_timestamps = 14

    def __init__(self):
        """
        Declares the variables used in solving the Heesch problem. The variables
        self.k_cor, self.shape, self.shape_size, self.shape_rad, self.grid and
        self.rotation_matrices should be defined in the initializer for the
        subclasses. Runs in O(1) time.
        """

        self.k_cor = None
        self.grid = None
        self.transforms = {}
        self.cells = {}
        self.rotation_matrices = []
        self.shape = None
        self.shape_size = None
        self.shape_rad = None
        self.sat = None
        self.model = None
        self.num_clauses = 0
        self.times = [0.0] * Heesch.num_timestamps
        self.ops = [0] * Heesch.num_timestamps

    @abstractmethod
    def generate_variables(self, start_corona: int):
        pass

    def construct_sat(self) -> Solver:
        # Constructs the SAT formula for the problem
        """
        Constructs the clauses of the boolean satisfiability expression and
        adds it to the SAT solver, self.sat. Generates the clauses as specified
        in the paper "Heesch Numbers of Unmarked Polyforms" by C. Kaplan. Uses
        a 'g3' solver in the library PySAT. Uses multiprocessing to generate
        certain types of clauses. Runs in O(mt^2) time, where t is the number of
        transformations and 2*m is the size of the shape. Note that t is bounded
        above by krn^2, so time complexity is equivalently O(k^2r^2mn^4).
        """

        # Note that a & b => c | d === -a | -b | c | d

        s = Solver(name="g3")

        self.times[3] = time()

        # 0-corona always used
        # O(1) time
        k_0 = self.transforms[
            (0, int(self.grid.size[0] / 2), int(self.grid.size[1] / 2), 0)
        ]
        s.add_clause([k_0[0]])
        self.num_clauses += 1

        self.times[4] = time()

        # if a transform is used, its cells are used
        # O(mt) time
        for k, v in self.transforms.items():
            # add each cell in the transform
            # O(m) time
            for c in v[1]:
                s.add_clause([-v[0], self.cells[tuple(c)]])
                self.num_clauses += 1

        self.times[5] = time()

        # if a cell is used, some transform must use it
        # O(mn^2t) time
        for k1, v1 in self.cells.items():
            lst = [-v1]

            # check if the cell is in the transform for
            # each transform
            # O(tm) time
            for k2, v2 in self.transforms.items():
                if list(k1) in v2[1].tolist():
                    lst.append(v2[0])
            s.add_clause(lst)
            self.num_clauses += 1

        self.times[6] = time()

        # if a transform is used in an interior corona, its halo cells must
        # be used
        # O(mt) time
        for k, v in self.transforms.items():
            # do not consider nth corona
            if k[0] == self.k_cor:
                continue

            # consider all cells in halo of transform
            # O(m) time
            for c in v[2]:
                s.add_clause([-v[0], self.cells[tuple(c)]])
                self.num_clauses += 1

        self.times[7] = time()

        # used transforms cannot overlap
        # O(mt^2/n) time
        with Pool(processes=cpu_count()) as pool:
            # takes O(mt^2/n) for each process
            # adding clause from list does not change time complexity
            for lst in pool.imap_unordered(
                self.check_overlap_mp, range(0, self.grid.size[0])
            ):
                for l in lst:
                    s.add_clause(l)
                    self.num_clauses += 1

        self.times[8] = time()

        # if a transform is used in a k corona, it must be adjacent to one
        # in a k-1 corona
        # O(mt^2/kn) time
        with Pool(processes=cpu_count()) as pool:
            # takes O(mt^2/kn) for each process
            # adding clause from list does not change time complexity
            for lst in pool.imap_unordered(
                self.check_halo_overlap_mp, range(0, self.grid.size[0])
            ):
                for l in lst:
                    s.add_clause(l)
                    self.num_clauses += 1

        self.times[9] = time()

        # if a transform is used in a k corona, it cannot be adjacent to one
        # in a 0...k-2 corona
        # O(mt^2) time
        for k1, v1 in self.transforms.items():
            for k2, v2 in self.transforms.items():
                # only consider coronas < k-1
                if k2[0] >= k1[0] - 1:
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                if self.grid.is_overlapping(v1[1], v2[2]):  # O(m)
                    s.add_clause([-v1[0], -v2[0]])
                    self.num_clauses += 1

        self.times[10] = time()

        # # TODO - add hole suppression
        # # not sure why we can't just have a and b and c and d => e
        # # where a, b, c, d, are the tiles adjacent to e
        # for k, v in self.cells.items():
        #     lst = []
        #
        #     # find all adjacent cells
        #     for c in self.grid.haloIdx(k):
        #         lst.append(-self.cells[tuple(c)])
        #
        #     lst.append(v)
        #     s.add_clause(lst)
        #     self.num_clauses += 1

        self.times[11] = time()

        self.sat = s
        return self.sat

    def solve_sat(self) -> list:
        """
        Solves the boolean satisfiability model given by self.sat. Uses the
        PySAT library. Return the model containing the variable specifications
        that make the boolean formula true. The integer values in the model
        correspond to the variable indices given in the self.transforms and
        self.cells dictionaries. Runs in non-polynomial time.
        """

        self.times[12] = time()

        if self.sat.solve():
            self.model = self.sat.get_model()
        else:
            self.model = None

        self.times[13] = time()

        return self.model

    def clear_sat(self):
        self.sat.delete()
        self.sat = None
        self.model = None

    def get_transform(self, idx: int):
        """
        Returns the transformation corresponding to the given index. The index
        given should correspond to the index found in the self.transforms
        dictionary. The key returned is a 4-tuple specified as (corona,
        x-translate, y-translate, rotation). Runs in O(1) time.
        """

        idx -= 1
        if idx >= (self.k_cor + 1) * self.grid.size[0] * self.grid.size[1] * len(
            self.rotation_matrices
        ):
            return None

        v1 = int(
            idx / (self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices))
        )
        idx -= v1 * (
            self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)
        )
        v2 = int(idx / (self.grid.size[1] * len(self.rotation_matrices)))
        idx -= v2 * (self.grid.size[1] * len(self.rotation_matrices))
        v3 = int(idx / (len(self.rotation_matrices)))
        idx -= v3 * (len(self.rotation_matrices))

        return v1, v2, v3, idx

    def get_transform_idx(self, key: tuple[int, int, int, int]):
        """
        Returns the index corresponding to the given key. The key should be a
        4-tuple specified as (corona, x-translate, y-translate, rotation). The
        index returned corresponds to the index found in the self.transforms
        dictionary. Runs in O(1) time.
        """

        return (
            key[0] * self.grid.size[0] * self.grid.size[1] * len(self.rotation_matrices)
            + key[1] * self.grid.size[1] * len(self.rotation_matrices)
            + key[2] * len(self.rotation_matrices)
            + key[3]
            + 1
        )

    def check_overlap_mp(self, i1: int):
        """
        Checks whether pairs of transforms overlap with each other. Only check
        pairs where one of the transforms is in the specified first grid index,
        since other processes will evaluate for other grid indices. Is meant
        to be used as part of a pool of processes. Returns a list of clauses in
        the form: [
            [ x_1, -x_2 ],
            [ x_3,  x_4 ]
            ...
        ] where each x_1, ..., x_4 is a variable. Runs in O(mt^2/n) time, where
        t is the number of transformations, 2*m is the shape size, and n*n is
        the size of the grid.
        """

        out = []
        # O(mt^2/n) time
        for (k1, v1), (k2, v2) in itertools.combinations(self.transforms.items(), 2):
            # for k1, v1 in self.transforms.items():

            # for k2 in range(0, self.k_cor):
            #     for k5 in range(0, len(self.rotation_matrices)):
            #         for k4 in range(k1[2], 0, -1):
            #             key = (k2, i1, k4, k5)
            #             if key not in self.transforms:
            #                 continue

            # v2 = self.transforms.get(key)

            if k1[1] != i1 and k2[1] != i1:
                continue

            if v1[0] >= v2[0]:
                continue

            # O(1) time
            if (
                max(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1]))
                > 2 * self.shape_rad
            ):
                continue

            # O(1) time
            # if min(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * self.shape_rad:
            #     break

            # checks if any two rows in the transform are the same
            # O(m) time
            # print(len(v1), len(v2))
            if len(v1) != 5 or len(v2) != 5:
                print(v1, v2)
                raise Exception(f'akdslkajdlsakd \n {len(v1)}, {len(v2)}, {v1}, {v2}')
            if self.grid.is_overlapping_set(v1[3], v2[3]):  # O(m) time
                out.append([-v1[0], -v2[0]])

                # for k4 in range(k1[2], self.grid.size[1], 1):
                #     key = (k2, i1, k4, k5)
                #     if key not in self.transforms:
                #         continue
                #
                #     v2 = self.transforms.get(key)
                #
                #     # if k1[1] != i1 and k2[1] != i1:
                #     #     continue
                #
                #     # if v1[0] >= v2[0]:
                #     #     continue
                #
                #     # O(1) time
                #     if max(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * self.shape_rad:
                #         continue
                #
                #     # O(1) time
                #     if min(abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])) > 2 * self.shape_rad:
                #         break
                #
                #     # checks if any two rows in the transform are the same
                #     # O(m) time
                #     if self.grid.is_overlapping(v1[1], v2[1]):  # O(m) time
                #         out.append([-v1[0], -v2[0]])

        return out

    def check_halo_overlap_mp(self, i1: int):
        """
        Checks whether a tranformation and the halo of another transformation
        overlap, where one transformation is in the specified first grid index,
        since the other processes will evaluate for other grid indices. Is meant
        to be used as part of a pool of processes. Returns a list of clauses in
        the form: [
            [ x_1, -x_2, ... ],
            [ x_3,  x_4, ... ]
            ...
        ] where each x_1, ..., x_4 is a variable. Runs in O(mt^2/kn) time, where
        t is the number of transformations, 2*m is the shape size, n*n is
        the size of the grid, and k is the number of coronas
        """

        # O(t^2m/kn) time
        out = []
        for k1, v1 in self.transforms.items():
            if k1[1] != i1:
                continue
            if k1[0] == 0:
                continue

            # iterate through the other transforms
            # O(tm/k) time
            l = [-v1[0]]
            for k2, v2 in self.transforms.items():
                # only consider the k and k-1 coronas
                if k1[0] != k2[0] + 1:
                    continue

                # Do not consider to be overlapping if more than 2*shape_radius away from each other
                # O(1) time
                if max(
                    abs(v1[1][0][0] - v2[1][0][0]), abs(v1[1][0][1] - v2[1][0][1])
                ) > 2 * (self.shape_rad + 1):
                    continue

                # checks if any two rows in the transform and second transform's
                # halo are the same
                # O(m) time
                if self.grid.is_overlapping_set(v1[3], v2[4]):  # O(m) time
                    l.append(v2[0])

            if len(l) > 0:
                out.append(l)

        return out

    def write(self, directory=".", plot=True):
        """
        Writes the output of the Heesch problem to a file. The filename is
        given by the time of writing. The directory parameter specifies the
        location the file should be written to, by default the current
        directory. The plot argument specifies whether to create a
        matplotlib figure of the corresponding tiling that solves the Heesch
        problem. Plotting requires the abstract function plot be implemented
        in a subclass.
        """

        filename = str(int(time()))

        with open(directory + "/" + filename + ".txt", "w+") as f:
            f.write(f"File: {filename}\n")
            f.write(f"Start Time: {self.times[0]}\n")
            f.write(f"Grid Size: {self.grid.size}\n")
            f.write(f"Coronas: {self.k_cor}\n")
            f.write(f"Shape Size: {self.shape_size}\n")
            f.write(f"Shape Radius: {self.shape_rad}\n")
            f.write(f"Shape: \n")
            f.write(f"{str(self.shape)}\n")

            f.write("--------------------------------------------------\n")

            f.write(f"Total Clauses: {self.num_clauses}\n")
            f.write(f"SAT Model: \n")
            if self.model is not None:
                f.write(f"\t{str(self.model)}\n")
            else:
                f.write(f"\tNo Model (Unsolvable)\n")
            f.write("Times: \n")
            sum_time = 0.0
            for i, t in enumerate(self.times):
                f.write(f'\tTimestamp {i:02d}: {t - sum_time - self.times[0]}s '
                        f'({self.ops[i] - self.ops[i - 1] if i != 0 else 0} ops)\n')

                sum_time = t - self.times[0]
            f.write(f"Total Time: {sum_time}s ({self.ops[-1]} ops)\n")

            f.write("--------------------------------------------------\n")

            f.write(f'Transforms Generated ({len(self.transforms)}): \n')
            k_count = [0] * (self.k_cor + 1)
            for key in self.transforms.keys():
                k_count[key[0]] += 1
            for k_idx, kc in enumerate(k_count):
                f.write(f'\tCorona {k_idx}: {kc}\n')

            if self.model is not None:
                num_transforms = len([i for i in self.model if i > 0 and self.get_transform(i) is not None])
                f.write(f'Transforms Used ({num_transforms}): \n')

                for i in self.model:
                    if i <= 0:
                        continue
                    t = self.get_transform(i)
                    if t is None:
                        break
                    f.write(f"Transform ID: {i} {self.get_transform(i)}: \n")
                    f.write(f"{str(self.transforms[t][1])}\n")
            else:
                f.write("\tNo Model (Unsolvable)")

        if plot:
            self.plot(show=False, write=True, filename=filename, directory=directory)

    @abstractmethod
    def plot(self, show, write, filename, directory):
        pass
