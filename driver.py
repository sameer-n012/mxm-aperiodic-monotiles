import glob
import os
import re
import sys
from math import floor
import numpy as np
from HeeschSat.PolykiteHeesch import PolykiteHeesch


def read_from_file(path: str, start_line: int = 0, end_line: int = -1):
    # TODO: transform (a,b) to (x,y,k) using above equations
    # and then try to read in the file line by line
    # and maybe we can start reading by a specific line
    # and stop if we need to

    """
    The coordinate model is related to that for polyhexes.  Each hexagon
    in $(6^3)$, oriented as in the coordinate model for polyhexes, is
    divided into six kites; consider those to be numbered $0$ through~$5$
    counterclockwise, with numbers $0$ through~$2$ being in the bottom
    half of the hexagon.  Then, in hexagon $(x, y)$, kite number~$k$ is
    given polykite coordinates $(3x + (k mod 3), 2y + floor(k/3))$.
    """

    kite_points = []
    lineCount = 0
    with open('results/all_kites.txt', "w+") as fout:

        for filename in glob.glob(os.path.join(path, "*.txt")):
            with open(os.path.join(os.getcwd(), filename), "r") as f:  # open in readonly mode
                # do your stuff
                # Using readlines()
                # startValue = 0
                # endValue = 2
                Lines = f.readlines()
                # Strips the newline character
                # print("New File")

                for line in Lines:
                    if lineCount < start_line:
                        lineCount += 1
                        continue
                    if end_line is not None and lineCount >= end_line:
                        return

                    dataList = [
                        tuple(map(int, match.group(1).split(",")))
                        for match in re.finditer(r"\((.*?)\)", line)
                    ]

                    converted_data = []
                    for a, b in dataList:

                        x = floor(a / 3)
                        y = floor(b / 2)
                        k = floor((b % 2) * 3 + (a % 3)) + 3
                        if k > 5:
                            k %= 3

                        converted_data.append([x, y, k])

                    lineCount += 1

                    cor = 1
                    while True:
                        ph = PolykiteHeesch(np.array(converted_data), coronas=cor)
                        print(f'{lineCount}. trying {cor}-corona for {converted_data}...')
                        ph.generate_variables()
                        ph.construct_sat()
                        ph.solve_sat()
                        if ph.model is None:
                            break
                        if ph.model is not None and cor >= 2:
                            ph.write(directory='tests/out', plot=True)
                        cor += 1
                    fout.write(f'({cor - 1}): {converted_data}\n')


if __name__ == '__main__':
    read_from_file(path="./data/", end_line=200)
