from HeeschSat.Heesch import Heesch
from HeeschSat.HexGrid import HexGrid
from numpy import identity
import itertools


class PolyhexagonHeesch(Heesch):

    def __init__(self, coronas):
        super().__init__()
        self.grid = HexGrid((10, 10))
        self.transforms = {}
        self.rotation_matrices = []
        self.k_cor = coronas

    def construct_sat(self):
        pass

    def solve_sat(self, sat: str):
        pass
