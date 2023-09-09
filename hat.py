from math import sqrt, pi, sin, cos

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

import numpy as np

class Tile:
    """
    A class representing tiles.
    
    """
    
    def __init__(self, vertices):
        """
        Initialize the tile with a numpy array of vertices.

        """

        self.vertices = vertices

    def draw(self):
        """
        Draws the tile.

        """

        fig, ax = plt.subplots()

        # Create Polygon object using this tile's vertices
        polygon = Polygon(self.vertices)
        patch = [polygon]

        p = PatchCollection(patch)
        ax.add_collection(p)

        plt.show()

        
class Hat(Tile):
    """
    A class representing the hat category of aperiodic monotiles.

    """

    def __init__(self, a, b, origin=[0, 0], reflected=False):
        """
        Initialize a reflected or non-reflected hat by creating vertices 
        according to the provided lengths. Hats are initialized with a fixed 
        orientation.  

        """

        # The angles of the hat's sides don't change, even with different 
        # side lengths.
        theta = pi / 3
        phi = pi / 6

        # Set the vertices of the hat as follows, according to the given origin
        # and side lengths.
        x0 = origin[0]
        y0 = origin[1]
        x1 = x0 + b
        y1 = y0
        x2 = x1 + b * cos(theta)
        y2 =  y1 + b * sin(theta)
        x3 = x2 + a * cos(phi)
        y3 = y2 - a * sin(phi)
        x4 = x3 + a * cos(phi)
        y4 = y2
        x5 = x4 + b * cos(theta)
        y5 = y0
        x6 = x4
        y6 = y5 - b * sin(theta)
        x7 = x6 + a * cos(phi)
        y7 = y6 - a * sin(phi)
        x8 = x7
        y8 = y7 - a
        x9 = x8 - b
        y9 = y8
        x10 = x9 - b * cos(theta)
        y10 = y9 + b * sin(theta)
        x11 = x10 - a * cos(phi)
        y11 = y10 - a * sin(phi)
        x12 = x0
        y12 = y0 - a

        vertices = np.array([[x0, y0], [x1, y1], [x2, y2], [x3, y3], [x4, y4],
                             [x5, y5], [x6, y6], [x7, y7], [x8, y8], [x9, y9],
                             [x10, y10], [x11, y11], [x12, y12]])

        if reflected:
            vertices = vertices.dot([[-1, 0], [0, 1]])

        # Calls the Tile constructor with the new 'vertices' array as input.
        super().__init__(vertices)


def draw_hat():
    hat = Hat(1, sqrt(3))
    hat.draw


def draw_reflected_hat():
    hat = Hat(1, sqrt(3), reflected=True)
    hat.draw()

