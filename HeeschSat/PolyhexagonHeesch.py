from HeeschSat.Heesch import Heesch
from HeeschSat.HexGrid import HexGrid
import numpy as np
from numpy import identity
import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon


class PolyhexagonHeesch(Heesch):

    def __init__(self, shape, coronas, grid_size):
        super().__init__(coronas)
        self.grid = HexGrid(grid_size)
        self.rotation_matrices = [np.matmul(np.linalg.matrix_power(
            np.array([[1, 1],
                      [0, -1]]), j),
            np.linalg.matrix_power(
                np.array([[0, -1],
                          [1, 1]]), i))
            for i in range(0, 6) for j in range(0, 2)
        ]
        self.k_cor = coronas
        self.shape = shape
        self.shape_size = self.shape.shape

    def plot(self):
        if self.model is None:
            return

        shapes = []
        for i in self.model:
            if i <= 0:
                continue
            t = self.get_transform(i)
            if t is None:
                break
            shapes.append(self.transforms[t][1])

        shapes = np.array(shapes)

        fig, ax = plt.subplots(1)
        ax.set_aspect("equal")

        i_s = self.grid.indices()
        i_ts = np.matmul(self.grid.basis, i_s.T).T
        for i, i_t in enumerate(i_ts):

            hexagon = RegularPolygon(
                (i_t[0], i_t[1]),
                numVertices=6,
                radius=np.sqrt(1 / 3),
                # orientation=np.pi / 6,
                alpha=1.0,
                facecolor='w',
                edgecolor="k",
                linewidth=0.2,
            )
            ax.add_patch(hexagon)
            ax.text(i_t[0], i_t[1], f'{i_s[i][0]}, {i_s[i][1]}',
                    verticalalignment='center',
                    horizontalalignment='center',
                    clip_on=True
            )

        for i, s in enumerate(shapes):
            s_t = np.matmul(self.grid.basis, s.T).T

            for j, c in enumerate(s_t):
                # alpha is transparency
                hexagon = RegularPolygon(
                    (c[0], c[1]),
                    numVertices=6,
                    radius=np.sqrt(1 / 3),
                    # orientation=np.pi / 6,
                    alpha=0.5,
                    facecolor=Heesch.plot_colors[i%len(Heesch.plot_colors)],
                    edgecolor="k",
                )
                ax.add_patch(hexagon)
                # ax.text(c[0], c[1], f'{s[j][0]}, {s[j][1]}', verticalalignment='center', '')
        plt.autoscale(enable=True)
        plt.show()


