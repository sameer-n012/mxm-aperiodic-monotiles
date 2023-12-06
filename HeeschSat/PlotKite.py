import os, glob
import numpy as np
import ast
import re
from math import floor

kite_points = []
path = "../JosephMyers/"
for filename in glob.glob(os.path.join(path, "*.txt")):
    with open(os.path.join(os.getcwd(), filename), "r") as f:  # open in readonly mode
        # do your stuff
        # Using readlines()
        startValue = 0
        endValue = 2
        Lines = f.readlines()
        lineCount = 0
        # Strips the newline character
        print("New File")
        for line in Lines:
            if lineCount not in range(startValue, endValue):
                lineCount += 1
                continue
            dataList = [
                tuple(map(int, match.group(1).split(",")))
                for match in re.finditer(r"\((.*?)\)", line)
            ]
            for a, b in dataList:
                x = floor(a / 3)
                y = floor(b / 2)
                k = floor((b % 2) * 3 + (a % 3)) + 3
                if k > 5:
                    k %= 3

                print(x, y, k)
            lineCount += 1


"""
The coordinate model is related to that for polyhexes.  Each hexagon
in $(6^3)$, oriented as in the coordinate model for polyhexes, is
divided into six kites; consider those to be numbered $0$ through~$5$
counterclockwise, with numbers $0$ through~$2$ being in the bottom
half of the hexagon.  Then, in hexagon $(x, y)$, kite number~$k$ is
given polykite coordinates $(3x + (k mod 3), 2y + floor(k/3))$.
    
"""


# TODO: transform (a,b) to (x,y,k) using above equations
# and then try to read in the file line by line
# and maybe we can start reading by a specific line
# and stop if we need to
