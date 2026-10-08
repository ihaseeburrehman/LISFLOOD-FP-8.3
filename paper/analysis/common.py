"""Shared helpers for results.py and figures.py (numpy and pandas only)."""
import glob
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
CASES = os.path.join(PAPER, "cases")
G = 9.81


def header(path):
    h = {}
    with open(path) as f:
        for _ in range(6):
            k, v = f.readline().split()
            h[k.lower()] = float(v)
    h["ncols"], h["nrows"] = int(h["ncols"]), int(h["nrows"])
    return h


def grid(path):
    """An ESRI ASCII grid as a float32 array (rows north to south)."""
    h = header(path)
    a = pd.read_csv(path, skiprows=6, sep=r"\s+", header=None, dtype=np.float32)
    return a.to_numpy(np.float32).reshape(h["nrows"], h["ncols"])


def steps(run_dir):
    """Number of time steps, from the last line of the .mass file."""
    m = glob.glob(os.path.join(run_dir, "*.mass"))
    with open(m[0]) as f:
        return int(f.read().strip().splitlines()[-1].split()[3])


def loop_time(run_dir):
    """Run time of the time loop in seconds, as LISFLOOD-FP reports it."""
    with open(os.path.join(run_dir, "run.log")) as f:
        t = [l.split()[2] for l in f if l.startswith("loop time")]
    return float(t[-1])


def outputs(run_dir, ext):
    """Output grids <resroot>-NNNN.<ext>, by output number."""
    out = {}
    for p in glob.glob(os.path.join(run_dir, "*-[0-9][0-9][0-9][0-9]." + ext)):
        out[int(re.search(r"-(\d{4})\." + ext + "$", p).group(1))] = p
    return out


def ritter(x, t, x0=1000.0, h0=1.0):
    """Exact depth of the dam break on a dry, frictionless bed (Ritter 1892)."""
    c0 = np.sqrt(G * h0)
    xi = x - x0
    h = np.where(xi <= -c0 * t, h0, (2 * c0 - xi / max(t, 1e-9)) ** 2 / (9 * G))
    return np.where(xi >= 2 * c0 * t, 0.0, h)


# Alzette river gauges: EPSG:2169 position, width of the surveyed section (m),
# and the downstream unit vector of the river there (taken once from the river
# centreline, as in the companion studies).
GAUGES = {
    "Walferdange": (77256, 81571, 30.0, 0.01176865738325593, 0.9999307469537055),
    "Steinsel": (77432, 82659, 50.0, 0.4195660314293628, 0.9077248180316626),
    "Pfaffenthal": (77409, 76226, 30.0, 0.15308181225518291, 0.9882135188089005),
    "Livange": (76151, 65753, 25.0, 0.9722263381555832, 0.23404261876116805),
    "Hesperange": (78623, 72404, 25.0, 0.7691691193433822, 0.6390452768376634),
}


def _transect_cells(h, x, y, W, px, py):
    """Cells a transect of width W crosses, with the exact length in each."""
    c = h["cellsize"]
    top = h["yllcorner"] + h["nrows"] * c
    t0, t1 = -W / 2.0, W / 2.0
    breaks = {t0, t1}
    if abs(px) > 1e-12:
        ka = int(np.floor(min(x + t0 * px, x + t1 * px) / c - h["xllcorner"] / c))
        kb = int(np.ceil(max(x + t0 * px, x + t1 * px) / c - h["xllcorner"] / c))
        for k in range(ka, kb + 1):
            t = (h["xllcorner"] + k * c - x) / px
            if t0 < t < t1:
                breaks.add(float(t))
    if abs(py) > 1e-12:
        ka = int(np.floor((top - max(y + t0 * py, y + t1 * py)) / c))
        kb = int(np.ceil((top - min(y + t0 * py, y + t1 * py)) / c))
        for k in range(ka, kb + 1):
            t = (top - k * c - y) / py
            if t0 < t < t1:
                breaks.add(float(t))
    ts = sorted(breaks)
    out = []
    for ta, tb in zip(ts[:-1], ts[1:]):
        if tb - ta <= 1e-12:
            continue
        tm = 0.5 * (ta + tb)
        out.append((int(np.floor((top - (y + tm * py)) / c)),
                    int(np.floor((x + tm * px - h["xllcorner"]) / c)), tb - ta))
    return out


def gauges(run_dir):
    """Discharge (m3/s) and depth (m) at the five gauges for every hourly output.

    Discharge is the integral of the unit discharge (Qx, Qy) normal to the
    surveyed section across the river; depth is the value in the gauge cell.
    """
    qx_files = outputs(run_dir, "Qx")
    h = header(next(iter(qx_files.values()))[:-3] + ".wd")
    top = h["yllcorner"] + h["nrows"] * h["cellsize"]
    rows = []
    for k in sorted(qx_files):
        QX, QY = grid(qx_files[k]), grid(qx_files[k][:-3] + ".Qy")
        WD = grid(qx_files[k][:-3] + ".wd")
        rec = {"output": k}
        for name, (x, y, W, ux, uy) in GAUGES.items():
            Q = sum((float(QX[r, c]) * ux + float(QY[r, c]) * uy) * L
                    for r, c, L in _transect_cells(h, x, y, W, -uy, ux))
            rec[name + "_Q"] = Q
            rec[name + "_depth"] = float(WD[int((top - y) // h["cellsize"]),
                                            int((x - h["xllcorner"]) // h["cellsize"])])
        rows.append(rec)
    return pd.DataFrame(rows)
