import matplotlib.pyplot as plt


def plot_kite(points):
    # Unpack points
    x = [point[0] for point in points]
    y = [point[1] for point in points]

    # Plot the kite shape
    plt.plot(x + [x[0]], y + [y[0]], marker="o")

    # Connect the points to form the kite
    plt.plot([x[0], x[2]], [y[0], y[2]], linestyle="-", color="black")
    plt.plot([x[1], x[3]], [y[1], y[3]], linestyle="-", color="black")

    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.title("Kite Plot")
    plt.grid(True)
    plt.show()


# Data points for the kite
kite_points = [(0, 0), (1, 0), (2, 1), (2, 0)]

"""
    
    The coordinate model is related to that for polyhexes.  Each hexagon
in $(6^3)$, oriented as in the coordinate model for polyhexes, is
divided into six kites; consider those to be numbered $0$ through~$5$
counterclockwise, with numbers $0$ through~$2$ being in the bottom
half of the hexagon.  Then, in hexagon $(x, y)$, kite number~$k$ is
given polykite coordinates $(3x + (k mod 3), 2y + floor(k/3))$.
    
"""
a = 0
b = 1
 # p = [a, b]

x = a / 3
y = b / 2
k = (b % 2) * 3 + (a % 3)


# TODO: transform (a,b) to (x,y,k) using above equations
# and then try to read in the file line by line
# and maybe we can start reading by a specific line
# and stop if we need to

[[x1, y1, k1], [x2, y2, k2], ...]  # numpy array

# Plot the kite
plot_kite(kite_points)
