from HeeschSat.Heesch import Heesch
from HeeschSat.HexGrid import HexGrid
import numpy as np
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
        self.shape_rad = max(self.shape_size[0], self.shape_size[1])

    def plot(self, show=True, write=False, filename=None, directory=None):
        if self.model is None:
            return

        fig, ax = plt.subplots(1)
        ax.set_aspect("equal")

        bounds_t = self.grid.apply_basis(np.array([list(self.grid.size)])).reshape(1, 2)[0]

        i_s = np.array([(i, j) for j in range(-self.grid.size[1], 2*self.grid.size[1])
                        for i in range(-self.grid.size[0], 2*self.grid.size[0])])
        i_ts = self.grid.apply_basis(i_s).T
        for i, i_t in enumerate(i_ts):
            if i_t[0] < 0*bounds_t[0] or i_t[1] < 0*bounds_t[0] or \
                    i_t[0] > 1*bounds_t[0] or i_t[1] > 1*bounds_t[1]:
                continue

            # grid hexagons
            hexagon = RegularPolygon(
                (i_t[0], i_t[1]),
                numVertices=6,
                radius=np.sqrt(1 / 3),
                # orientation=np.pi / 6,
                alpha=1.0,
                facecolor='w',
                edgecolor="k",
                # linewidth=None,
                linewidth=0.2,
                zorder=1.0
            )
            ax.add_patch(hexagon)

            # grid coords
            ax.text(i_t[0], i_t[1], f'{i_s[i][0]}, {i_s[i][1]}',
                    verticalalignment='center',
                    horizontalalignment='center',
                    clip_on=True,
                    fontsize=60 / self.grid.size[0],
                    zorder=5.0
                    )

        shapes = []
        trans = []
        colors = []
        c_idx = [0]*len(Heesch.c_majors)
        for i in self.model:
            if i <= 0:
                continue
            t = self.get_transform(i)
            if t is None:
                break

            # manage coloring of coronas
            colors.append(
                Heesch.plot_colors[
                    (Heesch.c_majors[
                        (Heesch.c_spacing*t[0]) % len(Heesch.c_majors)
                    ] + 2*c_idx[(Heesch.c_spacing*t[0]) % len(c_idx)])
                    % len(Heesch.plot_colors)
                ]
            )
            c_idx[(Heesch.c_spacing*t[0]) % len(c_idx)] += 1

            trans.append(t)
            shapes.append(self.transforms[t][1])

        shapes = np.array(shapes)

        # print(colors)

        for i, s in enumerate(shapes):
            s_t = self.grid.apply_basis(s).T
            # s_t = np.matmul(self.grid.basis, s.T).T

            for j, c in enumerate(s_t):

                # shape hexagons (filled, colored)
                # alpha is transparency
                # if trans[i][0] % 2 == 1:
                hexagon = RegularPolygon(
                    (c[0], c[1]),
                    numVertices=6,
                    radius=np.sqrt(1 / 3),
                    # orientation=np.pi / 6,
                    # alpha=0.5,
                    color=colors[i],
                    # facecolor=colors[i],
                    # edgecolor="k",
                    # linestyle='',
                    linewidth=None,
                    hatch='///' if trans[i][0] % 2 == 1 else '',
                    zorder=2.0
                )
                ax.add_patch(hexagon)
                # else:
                #     hexagon = RegularPolygon(
                #         (c[0], c[1]),
                #         numVertices=6,
                #         radius=np.sqrt(1 / 3),
                #         # orientation=np.pi / 6,
                #         alpha=0.5,
                #         facecolor=colors[i],
                #         edgecolor="k",
                #     )
                #     ax.add_patch(hexagon)

        ax.axis("off")
        # ax.set_xlim(-1, bounds_t[0]+1)
        # ax.set_ylim(-1, bounds_t[1]+1)
        plt.autoscale(enable=True)

        if write and filename is not None:
            if directory is not None:
                plt.savefig(directory + '/' + filename + '.png',
                            bbox_inches='tight',
                            dpi=100*self.grid.size[0]
                            )
            else:
                plt.savefig(filename + '.png',
                            bbox_inches='tight',
                            dpi=100*self.grid.size[0]
                            )

        if show:
            plt.show()

        return
