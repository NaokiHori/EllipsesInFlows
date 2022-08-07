import os
import sys
import numpy as np


def main(bef):
    ly = 1.
    offset = 0.
    aft = np.zeros(bef.shape)
    aft[0, 0] = bef[0, 0]
    aft[0, 1] = bef[0, 1]
    nitems = bef.shape[0]
    for n in range(nitems-1):
        x0 = bef[n  , 0]
        y0 = bef[n  , 1]
        x1 = bef[n+1, 0]
        y1 = bef[n+1, 1]
        if y1 < y0:
            offset += ly
        aft[n+1, 0] = x1
        aft[n+1, 1] = y1+offset
    return aft

if __name__ == "__main__":
    bfname = sys.argv[1]
    afname = sys.argv[2]
    data = np.loadtxt(bfname, usecols=[1, 2])
    data = main(data)
    np.savetxt(afname, data, fmt="% .7f", delimiter=" ")

