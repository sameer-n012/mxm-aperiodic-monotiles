import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np

# list contains center of a hexagon
offCoord = [[-2, 0], [-1, 0], [0, 0], [0, 1], [0, 2]]

fig, ax = plt.subplots(1)
# set aspect ration equal
ax.set_aspect("equal")

for c in offCoord:
    # alpha is transparency
    hexagon = RegularPolygon(
        (c[0], c[1]),
        numVertices=6,
        radius=np.sqrt(1 / 3),
        # orientation=np.pi / 6,
        alpha=0.2,
        edgecolor="k",
    )
    ax.add_patch(hexagon)
plt.autoscale(enable=True)
plt.show()
