import os, glob

kite_points = []
path = "../JosephMyers/"
for filename in glob.glob(os.path.join(path, "*.txt")):
    with open(os.path.join(os.getcwd(), filename), "r") as f:  # open in readonly mode
        # do your stuff
        # Using readlines()
        Lines = f.readlines()

        count = 0
        # Strips the newline character
        print("New File")
        for line in Lines:
            count += 1
            t = line.strip()
            # Data points
            # kite_points.append(t.spilt(") "))
            # print(kite_points)
            print("Line{}: {}".format(count, line.strip()))


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
