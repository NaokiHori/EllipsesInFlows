import numpy as np
from matplotlib import pyplot
from matplotlib import patches


def shift_and_rotate(disp, angle, p0):
    # shift a point p0, then rotate it
    # return the resulting point p1
    p1 = dict()
    x = p0["x"]+disp["x"]
    y = p0["y"]+disp["y"]
    p1["x"] = np.cos(angle)*x - np.sin(angle)*y
    p1["y"] = np.sin(angle)*x + np.cos(angle)*y
    return p1


def rotate_and_shift(disp, angle, p0):
    # rotate a point (x0, y0), then shift it
    # return the resulting point p1
    p1 = dict()
    x = p0["x"]
    y = p0["y"]
    p1["x"] = disp["x"] + np.cos(angle)*x - np.sin(angle)*y
    p1["y"] = disp["y"] + np.sin(angle)*x + np.cos(angle)*y
    return p1


def compute_evolute(e, t):
    # compute center of evolute and radius of fitted circle
    evolute = dict()
    a = e["a"]
    b = e["b"]
    evolute["x"] = a*(1.-np.power(b/a, 2.))*np.power(np.cos(t), 3.)
    evolute["y"] = b*(1.-np.power(a/b, 2.))*np.power(np.sin(t), 3.)
    num = np.power(np.power(a*np.sin(t), 2.)+np.power(b*np.cos(t), 2.), 1.5)
    den = a*b
    evolute["r"] = num/den
    return evolute


def find_normal_t(e, p):
    # compute parameter t, with which a point (a*cos(t), b*sin(t))
    # and the target point p gives a normal vector to the given ellipse
    # ref: https://blog.chatfield.io/simple-method-for-distance-to-ellipse/
    t = 0.25*np.pi
    while True:
        xe = e["a"]*np.cos(t)
        ye = e["b"]*np.sin(t)
        evolute = compute_evolute(e, t)
        dxe = xe-evolute["x"]
        dye = ye-evolute["y"]
        dxp = abs(p["x"])-evolute["x"]
        dyp = abs(p["y"])-evolute["y"]
        norme = np.hypot(dxe, dye)
        normp = np.hypot(dxp, dyp)
        dc = norme*np.arcsin((dxe*dyp-dye*dxp)/(norme*normp))
        dt = dc/np.sqrt(e["a"]**2.+e["b"]**2.-xe**2.-ye**2.)
        t += dt
        # limit in the first quadrant
        t = min(t, 0.5*np.pi)
        t = max(t,        0.)
        if abs(dt) < 1.e-8:
            break
    # apply result to the other quadrants
    if p["x"] < 0.:
        t = -t+np.pi
    if p["y"] < 0.:
        t = -t
    return t


def fit_circles(e0, e1):
    # core function, fitting two circles
    # read the documentation
    for n in range(10):
        if n == 0:
            # initialise temporary evolute
            p0 = {"x": e0["x"], "y": e0["y"]}
            p1 = {"x": e1["x"], "y": e1["y"]}
        # forward coordinate transform
        disp = {"x": -e0["x"], "y": -e0["y"]}
        angle = -e0["theta"]
        p1_ = shift_and_rotate(disp, angle, p1)
        # forward coordinate transform
        disp = {"x": -e1["x"], "y": -e1["y"]}
        angle = -e1["theta"]
        p0_ = shift_and_rotate(disp, angle, p0)
        # find optimum parameter t for each ellipse
        e0_t = find_normal_t(e0, p1_)
        e1_t = find_normal_t(e1, p0_)
        # compute evolutes and their radii
        evolute0_ = compute_evolute(e0, e0_t)
        evolute1_ = compute_evolute(e1, e1_t)
        # backward coordinate transform
        disp = {"x": +e0["x"], "y": +e0["y"]}
        angle = +e0["theta"]
        p0 = rotate_and_shift(disp, angle, evolute0_)
        # backward coordinate transform
        disp = {"x": +e1["x"], "y": +e1["y"]}
        angle = +e1["theta"]
        p1 = rotate_and_shift(disp, angle, evolute1_)
    evolute0 = {"x": p0["x"], "y": p0["y"], "r": evolute0_["r"]}
    evolute1 = {"x": p1["x"], "y": p1["y"], "r": evolute1_["r"]}
    return evolute0, evolute1


def plot_result(e0, e1, c0, c1):
    # wrapper of matplotlib.patches to visualise
    # resulting ellipses and fitted circles
    def get_patch_of_ellipse(e):
        x = e["x"]
        y = e["y"]
        width = 2.*e["a"]
        height = 2.*e["b"]
        angle = 180./np.pi*e["theta"]
        keywords = {
                "xy": (x, y),
                "width": width,
                "height": height,
                "angle": angle,
                "fc": "none",
                "ec": "#FF0000",
        }
        return patches.Ellipse(**keywords)

    def get_patch_of_circle(c):
        x = c["x"]
        y = c["y"]
        radius = c["r"]
        keywords = {
                "xy": (x, y),
                "radius": radius,
                "fc": "none",
                "ec": "#0000FF",
        }
        return patches.Circle(**keywords)

    def compute_bounds(es, cs):
        # determine plot ranges
        # i.e., arguments of set_[xy]lim
        xmin = +np.inf
        xmax = -np.inf
        ymin = +np.inf
        ymax = -np.inf
        for e in es:
            xmin = min(xmin, e["x"]-max(e["a"], e["b"]))
            xmax = max(xmax, e["x"]+max(e["a"], e["b"]))
            ymin = min(ymin, e["y"]-max(e["a"], e["b"]))
            ymax = max(ymax, e["y"]+max(e["a"], e["b"]))
        for c in cs:
            xmin = min(xmin, c["x"]-c["r"])
            xmax = max(xmax, c["x"]+c["r"])
            ymin = min(ymin, c["y"]-c["r"])
            ymax = max(ymax, c["y"]+c["r"])
        margin = 0.25
        xmin -= margin
        xmax += margin
        ymin -= margin
        ymax += margin
        return (xmin, xmax), (ymin, ymax)
    # prepare matplotlib objects
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    # draw ellipses and circles
    ax.add_patch(get_patch_of_ellipse(e0))
    ax.add_patch(get_patch_of_ellipse(e1))
    ax.add_patch(get_patch_of_circle(c0))
    ax.add_patch(get_patch_of_circle(c1))
    # set auxiliary parameters
    xlim, ylim = compute_bounds([e0, e1], [c0, c1])
    keywords = {
            "xlim": xlim,
            "ylim": ylim,
            "aspect": "equal",
    }
    ax.set(**keywords)
    # visualise
    pyplot.show()
    # clean-up
    pyplot.close()


if __name__ == "__main__":
    # initialise ellipses
    e0 = {"a": 2.0, "b": 1.5, "x": 0.5, "y": 0.5, "theta": 0.2}
    e1 = {"a": 1.5, "b": 1.0, "x": 2.0, "y": 2.5, "theta": 2.0}
    # fit circles
    c0, c1 = fit_circles(e0, e1)
    # plot
    plot_result(e0, e1, c0, c1)
