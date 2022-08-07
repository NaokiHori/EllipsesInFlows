import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches as patches
from tqdm import tqdm


def load(dname):
    time = np.load("{}/time.npy".format(dname))
    lx   = np.load("{}/lx.npy".format(dname))
    ly   = np.load("{}/ly.npy".format(dname))
    Re   = np.load("{}/Re.npy".format(dname))
    xc   = np.load("{}/xc.npy".format(dname))
    yc   = np.load("{}/yc.npy".format(dname))
    ux   = np.load("{}/ux.npy".format(dname))
    uy   = np.load("{}/uy.npy".format(dname))
    pas  = np.load("{}/particle_as.npy".format(dname))
    pbs  = np.load("{}/particle_bs.npy".format(dname))
    pxs  = np.load("{}/particle_xs.npy".format(dname))
    pys  = np.load("{}/particle_ys.npy".format(dname))
    pazs = np.load("{}/particle_azs.npy".format(dname))
    return time, lx, ly, Re, xc, yc, ux, uy, pas, pbs, pxs, pys, pazs


if __name__ == "__main__":
    if len(sys.argv) == 2:
        movie = True
    else:
        movie = False
    root = "output/save"
    dnames = sorted(["{}/{}".format(root, dname) for dname in os.listdir(root)])
    for cnt, dname in enumerate(tqdm(dnames)):
        time, lx, ly, Re, x, y, ux, uy, pas, pbs, pxs, pys, pazs = load(dname)
        x = x[1:-1]
        ux = 0.5*(ux[:, 1:]+ux[:, :-1])
        uy = uy[:, 1:-1]
        uy = np.vstack([uy, uy[0, :]])
        uy = 0.5*(uy[1:, :]+uy[:-1, :])
        # x, y = np.meshgrid(x, y)
        vel = np.sqrt(np.power(ux, 2.)+np.power(uy, 2.))
        # visualise transposed array since my screen is wider
        fig = plt.figure(figsize=(9, 9))
        ax111 = fig.add_subplot(111)
        ax111.contourf(y, x, vel.T, cmap="gnuplot", levels=51)
        for pindex, (pa, pb, px, py, paz) in enumerate(zip(pas, pbs, pxs, pys, pazs)):
            for yshift in [-ly, 0., ly]:
                e = patches.Ellipse(xy=(py+yshift, px), width=2.*pa, height=2.*pb, angle=90.-180./np.pi*paz, fc="#AAAAAA", ec="#000000")
                ax111.add_patch(e)
        kwrds = {
                "title": "{:.2f} ({: .2f} - {: .2f})".format(time, np.min(vel), np.max(vel)),
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
        plt.savefig("tmp/images/image{:010d}.png".format(cnt), bbox_inches="tight", pad_inches=0, dpi=300)
        # if movie:
        #     plt.show(block=False)
        #     plt.pause(0.01)
        # else:
        #     plt.show()
        plt.close()

