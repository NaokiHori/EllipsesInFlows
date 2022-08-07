import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import patches as pts


def load(dname):
    lx   = np.load("{}/lx.npy".format(dname))
    ly   = np.load("{}/ly.npy".format(dname))
    xc   = np.load("{}/xc.npy".format(dname))
    yc   = np.load("{}/yc.npy".format(dname))
    ux   = np.load("{}/ux.npy".format(dname))
    uy   = np.load("{}/uy.npy".format(dname))
    pa   = np.load("{}/particle_as.npy".format(dname))[0]
    pb   = np.load("{}/particle_bs.npy".format(dname))[0]
    px   = np.load("{}/particle_xs.npy".format(dname))[0]
    py   = np.load("{}/particle_ys.npy".format(dname))[0]
    paz  = np.load("{}/particle_azs.npy".format(dname))[0]
    return lx, ly, xc, yc, ux, uy, pa, pb, px, py, paz


if __name__ == "__main__":
    root = sys.argv[1]
    # find last directory
    dnames = sorted(["{}/{}".format(root, dname) for dname in os.listdir(root)])
    dname = dnames[-1]
    # load all data
    lx, ly, x, y, ux, uy, pa, pb, px, py, paz = load(dname)
    # cut boundaries, interpolate to cell centers
    x = x[1:-1]
    ux = 0.5*(ux[:, 1:]+ux[:, :-1])
    uy = uy[:, 1:-1]
    uy = np.vstack([uy, uy[0, :]])
    uy = 0.5*(uy[1:, :]+uy[:-1, :])
    fig = plt.figure(figsize=(8, 8))
    ax111 = fig.add_subplot(111)
    ax111.contourf(y, x, uy.T, cmap="gnuplot", vmin=np.min(uy), vmax=np.max(uy), levels=51)
    for yshift in [-ly, 0., ly]:
        e = pts.Ellipse(xy=(py+yshift, px), width=2.*pa, height=2.*pb, angle=90.-180./np.pi*paz, fc="#AAAAAA", ec="#000000")
        ax111.add_patch(e)
    kwrds = {
            "title": "",
            "aspect": "equal",
            "xlim": [0., ly],
            "ylim": [0., lx],
            "xlabel": "",
            "ylabel": "",
            "xticks": [],
            "yticks": [],
    }
    ax111.set(**kwrds)
    ax111.spines["bottom"].set_visible(False)
    ax111.spines["top"].set_visible(False)
    ax111.spines["left"].set_visible(False)
    ax111.spines["right"].set_visible(False)
    plt.savefig("snapshot.png", bbox_inches="tight", pad_inches=0, dpi=150)
    plt.close()

