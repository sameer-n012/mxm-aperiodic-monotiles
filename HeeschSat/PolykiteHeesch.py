import math

from HeeschSat.Heesch import Heesch
from HeeschSat.KiteGrid import KiteGrid
import numpy as np
from time import time
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon


class PolykiteHeesch(Heesch):
    def __init__(self, shape, coronas=0, grid_size=None):
        super().__init__()
        self.k_cor = coronas
        self.shape = shape
        self.shape_size = self.shape.shape
        self.shape_rad = np.max(np.ptp(self.shape[:, 0:2], axis=0)) + 1
        if grid_size is None:
            grid_size = int(2 * self.shape_rad * (self.k_cor + 1))
            grid_size = (grid_size, grid_size, 6)
        self.grid = KiteGrid(grid_size, 0)
        self.num_rotations = None

    def generate_variables(self, start_corona: int = 0):
        # TODO implement, not changed much yet
        """
        Generates the two sets of variables that are used to solve the Heesch
        problem. These are the transformation variables (one variable for each
        unique possible transformation in each corona) and the cell variables
        (one variable for each cell in the grid). The 0-corona only contains a
        single transformation in the center of the grid. All other coronas
        contain transformations that are bounded by the shape radius and the
        corona number. Transformations outside the boundary of the grid are not
        generated. Runs in O(krm^2n^4) time, where k is the number of coronas, r
        is the number of rotation/reflection matrices, 2*m is the size of the
        shape, and n*n is the size of the grid.
        """

        self.times[0] = time()
        # midpoint of the grid
        mid = int(self.grid.size[0] / 2), int(self.grid.size[1] / 2)

        # get the used rotation matrices (some might result in rotational symmetry)
        # O(rm) time
        t_rotations = self.grid.get_transformations(self.shape)
        self.num_rotations = len(t_rotations)
        rotate_indices = self.check_rotational_symmetry()

        # calculate the maximum number of transform
        max_transforms = (
                (self.k_cor + 1)
                * self.grid.size[0]
                * self.grid.size[1]
                * self.num_rotations
        )

        # include single 0-corona transform as original shape
        # O(m^2n^2) time
        translate_mat = np.empty(self.shape_size, dtype=int)  # O(m) time
        translate_mat[:, 0] = mid[0]  # O(m) time
        translate_mat[:, 1] = mid[1]  # O(m) time
        translate_mat[:, 2] = 0
        transform = self.shape + translate_mat  # O(m) time
        halo = np.unique([i for c in transform for i in self.grid.haloIdx(c)], axis=0)  # O(m^2n^2) time
        transform_set = set(map(tuple, transform))
        halo_set = set(map(tuple, halo)) - transform_set

        key = (0, mid[0], mid[1], 0)
        print(key)
        print(self.num_rotations)
        print(self.get_transform_idx(key))
        self.transforms[key] = self.get_transform_idx(key), transform, halo, transform_set, halo_set

        corona_halos = [set() for _ in range(0, self.k_cor + 1)]
        corona_halos[0] = halo_set

        # transform variables (if a transform is used)
        # O(kn^4rm^2) time
        offset = self.grid.size[0] * self.grid.size[1] * self.num_rotations
        for idx, val in enumerate(itertools.product(range(start_corona, self.k_cor + 1),
                                                    range(0, self.grid.size[0]),
                                                    range(0, self.grid.size[1]),
                                                    range(0, len(t_rotations))
                                                    )):

            # do not generate any 0-corona transforms
            if val[0] == 0:
                continue

            if tuple(val) in self.transforms and self.transforms[tuple(val)] is not None:
                continue

            # do not generate any rotation transforms where rotation is
            # symmetric to another rotation that was previously generated
            if rotate_indices[val[3]] == 0:
                continue

            # build the transform by translating and rotating it
            # O(m) time
            translate_mat = np.empty(self.shape_size, dtype=int)  # O(m) time
            translate_mat[:, 0] = val[1]  # O(m) time
            translate_mat[:, 1] = val[2]  # O(m) time
            transform = t_rotations[val[3]] + translate_mat  # O(m) time

            # ignore any transforms that are more than k_cor*shape_radius from the center
            # O(m) time
            if max(abs(val[1] - mid[0]), abs(val[2] - mid[1])) > (val[0] + 1) * (
                    self.shape_rad + 0
            ):
                continue

            # ignore any transforms that are less than k_cor from the center
            # O(m) time
            if max(abs(val[1] - mid[0]), abs(val[2] - mid[1])) < (val[0] - 1):
                continue

            # only include transforms that fall within the bounds of the grid
            if not self.grid.in_bounds(transform):
                continue

            # adjacent = False
            # for i2, v2 in enumerate(itertools.product(range(0, self.grid.size[0]),
            #                                           range(0, self.grid.size[1]),
            #                                           range(0, len(self.rotation_matrices))
            #                                           )):
            #     if (val[0] - 1, v2[0], v2[1], v2[2]) in self.transforms:
            #         if self.grid.is_overlapping(self.transforms[(val[0] - 1, v2[0], v2[1], v2[2])][2], transform):
            #             adjacent = True
            #
            # if not adjacent:
            #     continue

            # maybe could generate halos of an entire corona
            # precompute set reps of each transform
            transform_set = set(map(tuple, transform))
            if not self.grid.is_overlapping_set(corona_halos[val[0] - 1], transform_set):
                continue

            # O(m^2n^2) time
            halo = np.unique([i for c in transform for i in self.grid.haloIdx(c)], axis=0)  # O(m^2n^2)
            halo_set = set(map(tuple, halo)) - transform_set
            corona_halos[val[0]] = corona_halos[val[0]] | halo_set
            self.transforms[tuple(val)] = (idx + offset + 1, transform, halo, transform_set, halo_set)

        self.times[1] = time()

        # cell variables (if a cell is taken)
        offset = max_transforms + 1
        for idx, i in enumerate(self.grid.indices()):
            self.cells[tuple(i)] = offset + idx

        self.times[2] = time()

    def check_rotational_symmetry(self):
        """
        Checks whether any rotational transformations to the given shape result
        in duplicate shapes. Returns a list with 0's in the indices of
        rotational transforms that contain duplicates and 1's in the indices of
        rotational transforms that are not duplicates. The indices correspond
        to the indices of self.rotation_matrices. Runs in O(rm^2) time where
        r is the number of rotation/reflection matrices and m*2 is the size
        of the shape.
        """

        seen = []
        t_rotations = self.grid.get_transformations(self.shape)
        out = [1] * len(t_rotations)

        for idx, rt in enumerate(t_rotations):
            for x in range(-self.shape_rad, self.shape_rad):
                for y in range(-self.shape_rad, self.shape_rad):
                    translate_mat = np.empty(self.shape_size, dtype=int)
                    translate_mat[:, 0] = x
                    translate_mat[:, 1] = y
                    translate_mat[:, 2] = 0
                    transform = rt + translate_mat

                    t = set(map(tuple, transform))
                    if t in seen:
                        out[idx] = 0

            if out[idx] == 1:
                seen.append(set(map(tuple, rt)))

        return out

    def get_transform(self, idx: int):
        """
        Returns the transformation corresponding to the given index. The index
        given should correspond to the index found in the self.transforms
        dictionary. The key returned is a 4-tuple specified as (corona,
        x-translate, y-translate, rotation). Runs in O(1) time.
        """

        idx -= 1
        if idx >= (self.k_cor + 1) * self.grid.size[0] * self.grid.size[1] * self.num_rotations:
            return None

        v1 = int(idx / (self.grid.size[0] * self.grid.size[1] * self.num_rotations))
        idx -= v1 * (self.grid.size[0] * self.grid.size[1] * self.num_rotations)
        v2 = int(idx / (self.grid.size[1] * self.num_rotations))
        idx -= v2 * (self.grid.size[1] * self.num_rotations)
        v3 = int(idx / self.num_rotations)
        idx -= v3 * self.num_rotations

        return v1, v2, v3, idx

    def get_transform_idx(self, key: tuple[int, int, int, int]):
        """
        Returns the index corresponding to the given key. The key should be a
        4-tuple specified as (corona, x-translate, y-translate, rotation). The
        index returned corresponds to the index found in the self.transforms
        dictionary. Runs in O(1) time.
        """

        return (
                key[0] * self.grid.size[0] * self.grid.size[1] * self.num_rotations
                + key[1] * self.grid.size[1] * self.num_rotations
                + key[2] * self.num_rotations
                + key[3]
                + 1
        )

    def plot(self, show=True, write=False, filename=None, directory=None):
        # TODO edit
        if self.model is None:
            return

        fig, ax = plt.subplots(1)
        ax.set_aspect("equal")

        print(np.array([list(self.grid.size)]))

        bounds_t = self.grid.apply_basis(np.array([list(self.grid.size)])).reshape(
            1, 3
        )[0]

        print('here')

        i_s = np.array(
            [
                (i, j, k)
                for j in range(-self.grid.size[1], 2 * self.grid.size[1])
                for i in range(-self.grid.size[0], 2 * self.grid.size[0])
                for k in range(0, 6)
            ]
        )
        i_ts = self.grid.apply_basis(i_s).T
        for i, i_t in enumerate(i_ts):
            if (
                    i_t[0] < 0 * bounds_t[0]
                    or i_t[1] < 0 * bounds_t[0]
                    or i_t[0] > 1 * bounds_t[0]
                    or i_t[1] > 1 * bounds_t[1]
            ):
                continue

            # grid hexagons
            hexagon = RegularPolygon(
                (i_t[0], i_t[1]),
                numVertices=6,
                radius=np.sqrt(1 / 3),
                # orientation=np.pi / 6,
                alpha=1.0,
                facecolor="w",
                edgecolor="k",
                # linewidth=None,
                linewidth=0.2,
                zorder=1.0,
            )
            ax.add_patch(hexagon)

            # grid coords
            ax.text(
                i_t[0],
                i_t[1],
                f"{i_s[i][0]}, {i_s[i][1]}",
                verticalalignment="center",
                horizontalalignment="center",
                clip_on=True,
                fontsize=60 / self.grid.size[0],
                zorder=5.0,
            )

        shapes = []
        trans = []
        colors = []
        c_idx = [0] * len(Heesch.c_majors)
        for i in self.model:
            if i <= 0:
                continue
            t = self.get_transform(i)
            if t is None:
                break

            # manage coloring of coronas
            colors.append(
                Heesch.plot_colors[
                    (
                            Heesch.c_majors[
                                (Heesch.c_spacing * t[0]) % len(Heesch.c_majors)
                                ]
                            + 2 * c_idx[(Heesch.c_spacing * t[0]) % len(c_idx)]
                    )
                    % len(Heesch.plot_colors)
                    ]
            )
            c_idx[(Heesch.c_spacing * t[0]) % len(c_idx)] += 1

            trans.append(t)
            shapes.append(self.transforms[t][1])

        shapes = np.array(shapes)

        # print(colors)

        # shape hexagons
        for i, s in enumerate(shapes):
            s_t = self.grid.apply_basis(s).T
            # s_t = np.matmul(self.grid.basis, s.T).T

            for j, c in enumerate(s_t):
                # shape hexagons (filled, colored)
                # alpha is transparency
                # if trans[i][0] % 2 == 1:

                x = [
                    c[0],
                    1.0/2 * math.cos((c[2] % 6) * 2 * math.pi / 6) + c[0],
                    math.sqrt(1.0 / 3) * math.cos((c[2] % 6) * 2 * math.pi / 6 + math.pi / 6) + c[0],
                    1.0/2 * math.cos(((c[2] + 1) % 6) * 2 * math.pi / 6) + c[0]
                ]

                y = [
                    c[1],
                    1.0/2 * math.sin((c[2] % 6) * 2 * math.pi / 6) + c[1],
                    math.sqrt(1.0 / 3) * math.sin((c[2] % 6) * 2 * math.pi / 6 + math.pi / 6) + c[1],
                    1.0/2 * math.sin(((c[2] + 1) % 6) * 2 * math.pi / 6) + c[1]
                ]

                ax.fill(x, y,
                        # alpha=0.5,
                        color=colors[i],
                        # facecolor=colors[i],
                        # edgecolor="k",
                        # linestyle='',
                        linewidth=1 / self.grid.size[0],
                        hatch="///" if trans[i][0] % 2 == 1 else "",
                        zorder=2.0,
                        )

        ax.axis("off")
        # ax.set_xlim(-1, bounds_t[0]+1)
        # ax.set_ylim(-1, bounds_t[1]+1)
        plt.autoscale(enable=True)

        if write and filename is not None:
            if directory is not None:
                plt.savefig(
                    directory + "/" + filename + ".png",
                    bbox_inches="tight",
                    dpi=100 * self.grid.size[0],
                )
            else:
                plt.savefig(
                    filename + ".png", bbox_inches="tight", dpi=100 * self.grid.size[0]
                )

        if show:
            plt.show()

        return
