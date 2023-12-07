import numpy as np
import itertools
from HeeschSat.PolyhexagonHeesch import PolyhexagonHeesch
from HeeschSat.HexGrid import HexGrid
import matplotlib.pyplot as plt




x = [0, 1, 2]
y = [0, 3, 0]
fig, ax = plt.subplots(1)
ax.fill(x, y, alpha=1.0,
                facecolor="b",
                edgecolor="k",
                # linewidth=None,
                linewidth=0.2,
                zorder=1.0,)
plt.show()
