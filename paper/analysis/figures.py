#!/usr/bin/env python3
"""Draw the figures of the paper from the runs of run_cases.sh.

    python3 paper/analysis/figures.py [WORK]

Writes fig1_dambreak, fig2_ea4, fig3_alzette, fig4_gauges and fig5_speed (pdf
and png) into WORK (default paper/work); a figure whose runs are missing is
skipped. Run results.py first (it writes the gauge table used by fig4).

Colours: CPU green, CUDA vermillion, Metal blue, reference black. Where results
coincide every line stays visible: Metal and CUDA use the same dash pattern
shifted by half a period, CPU is dotted on top, the reference thin and solid.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
from matplotlib.colors import LogNorm, to_rgb
from matplotlib.lines import Line2D

from common import CASES, GAUGES, grid, header, loop_time, outputs, ritter, steps

WORK = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "work"))
CPU, CUDA, METAL, REF = "#009E73", "#D55E00", "#0072B2", "#000000"
STYLE = {"metal": dict(color=METAL, lw=1.6, ls=(0, (4, 4)), zorder=2, label="Metal"),
         "cuda83": dict(color=CUDA, lw=1.6, ls=(4, (4, 4)), zorder=2, label="CUDA 8.3"),
         "cpu": dict(color=CPU, lw=1.6, ls=(0, (1, 2.5)), zorder=3, label="CPU")}
S_REF = dict(color=REF, lw=0.7, ls="-", zorder=4)
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "legend.handlelength": 4.0, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})


def run(case, label):
    d = os.path.join(WORK, f"{case}_{label}")
    return d if os.path.isdir(d) and outputs(d, "wd") else None


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(WORK, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def nolabel(s):
    return {k: v for k, v in s.items() if k != "label"}


def fig1():
    runs = [(l, run("dambreak", l)) for l in STYLE if run("dambreak", l)]
    if not runs:
        return
    x = (np.arange(400) + 0.5) * 5.0
    fig, ax = plt.subplots(figsize=(4.8, 2.9), constrained_layout=True)
    for k, t in ((1, 20), (2, 40), (3, 60)):
        first = k == 1
        ax.plot(x, ritter(x, t), **S_REF, label="exact (Ritter, 1892)" if first else None)
        for lab, d in runs:
            ax.plot(x, grid(outputs(d, "wd")[k])[1], **(STYLE[lab] if first else nolabel(STYLE[lab])))
        xl = 1000 + 1.2 * np.sqrt(9.81) * t
        ax.annotate(f"{t} s", (xl, float(ritter(np.array(xl), t))), xytext=(5, 3),
                    textcoords="offset points", fontsize=7)
    ax.set(xlim=(700, 1420), ylim=(0, 1.05), xlabel="distance x (m)", ylabel="water depth (m)")
    ax.legend(frameon=False, loc="upper right")
    save(fig, "fig1_dambreak")


def fig2():
    ref = os.path.join(CASES, "ea_test4", "reference-5m")
    runs = [(l, run("ea4_5m", l)) for l in STYLE if run("ea4_5m", l)]
    if not runs or not os.path.isdir(ref):
        return
    row = 199                      # the grid row just north of y = 1000 m, through the breach
    x = (np.arange(200) + 0.5) * 5.0
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True,
                               gridspec_kw=dict(width_ratios=[1.6, 1]))
    for h in (1, 3, 5):
        first = h == 1
        a.plot(x, grid(outputs(ref, "wd")[h])[row], **S_REF,
               label="published (Shaw et al., 2021)" if first else None)
        for lab, d in runs:
            a.plot(x, grid(outputs(d, "wd")[h])[row], **(STYLE[lab] if first else nolabel(STYLE[lab])))
    for h, (xl, yl) in {1: (330, 0.045), 3: (790, 0.07), 5: (850, 0.152)}.items():
        a.text(xl, yl, f"{h} h", fontsize=7)
    a.set(xlim=(0, 1000), ylim=(0, None), xlabel="distance from the breach x (m)",
          ylabel="water depth (m)", title="(a) depth through the breach")
    a.legend(frameon=False, loc="upper right")
    r1 = grid(outputs(ref, "wd")[1])[row]
    for lab, d in runs:
        b.plot(x, 1000 * (grid(outputs(d, "wd")[1])[row] - r1), **nolabel(STYLE[lab]))
    b.axhline(0, color=REF, lw=0.6)
    b.set(xlim=(0, 200), xlabel="distance from the breach x (m)",
          ylabel="difference from published (mm)", title="(b) near the breach, 1 h")
    save(fig, "fig2_ea4")


def fig3(k=7):
    m, c = run("alzette", "metal"), run("alzette", "cpu")
    dem = os.path.join(CASES, "alzette", "Alzette_sub_basin_10m_bridge_burn.asc")
    if not (m and c and os.path.exists(dem)):
        return
    M, C, h = grid(outputs(m, "wd")[k]), grid(outputs(c, "wd")[k]), header(outputs(m, "wd")[k])
    basin = grid(dem) > -9000
    ext = (h["xllcorner"], h["xllcorner"] + h["ncols"] * h["cellsize"],
           h["yllcorner"], h["yllcorner"] + h["nrows"] * h["cellsize"])
    WHITE, BASIN = np.ones(3), np.full(3, 0.92)
    norm = LogNorm(0.01, 10)
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 4.6), constrained_layout=True, sharey=True)
    # each panel is one plain RGB image, coloured cell by cell here
    rgb = np.empty(M.shape + (3,)); rgb[:] = WHITE; rgb[basin] = BASIN
    w = basin & (M > 0.01)
    rgb[w] = cm.Blues(norm(M[w]))[:, :3]
    ax[0].imshow(rgb, extent=ext, interpolation="nearest")
    for n, (x, y, *_) in GAUGES.items():
        ax[0].plot(x, y, "o", ms=4.5, mfc="white", mec="k", mew=0.8)
        ax[0].annotate(n, (x, y), xytext=(5, 0), textcoords="offset points", fontsize=6.5, va="center")
    ax[0].set_title("(a) water depth, Metal, 15 July 01:00")
    ax[0].set_ylabel("Northing, EPSG:2169 (m)")
    fig.colorbar(cm.ScalarMappable(norm=norm, cmap="Blues"), ax=ax[0], shrink=0.5, label="water depth (m)")
    # (b) at the printed resolution: each pixel is B x B cells and shows the
    # largest difference among them, so no differing cell is lost
    wet = basin & ((M > 0.01) | (C > 0.01))
    d = np.abs(M - C) * 1000
    classes = [((d > 0) & (d <= 1.5), "#E69F00", "1 mm"),
               ((d > 1.5) & (d <= 10.5), "#CC79A7", "2 to 10 mm"),
               (d > 10.5, "#000000", "more than 10 mm")]
    B = 4
    ny, nx = (M.shape[0] // B) * B, (M.shape[1] // B) * B
    blk = lambda a, f: f(a[:ny, :nx].reshape(ny // B, B, nx // B, B), axis=(1, 3))
    cls = np.zeros(M.shape, np.int8); handles = []
    for i, (sel, col, lab) in enumerate(classes, start=1):
        s = wet & sel; cls[s] = i
        handles.append(Line2D([], [], ls="", marker="s", ms=5, color=col, label=f"{lab} ({int(s.sum()):,} cells)"))
    cls_b, basin_b, wet_b = blk(cls, np.max), blk(basin, np.any), blk(wet, np.any)
    dep_b = blk(np.where(wet, np.maximum(M, C), 0.0), np.max)
    rgb = np.empty(cls_b.shape + (3,)); rgb[:] = WHITE; rgb[basin_b] = BASIN
    g = 0.80 - 0.35 * np.clip(norm(np.where(wet_b, dep_b, 0.01)), 0, 1)
    rgb[wet_b] = np.stack([g, g, g], axis=-1)[wet_b]
    for i, (_, col, _) in enumerate(classes, start=1):
        rgb[cls_b == i] = to_rgb(col)
    ax[1].imshow(rgb, extent=(ext[0], ext[0] + nx * h["cellsize"], ext[3] - ny * h["cellsize"], ext[3]),
                 interpolation="nearest")
    ax[1].set_title(f"(b) cells where Metal and CPU differ\n{int((wet & (d > 0)).sum()):,} of {int(wet.sum()):,} wet cells")
    ax[1].legend(handles=handles, frameon=False, loc="upper center", fontsize=7, bbox_to_anchor=(0.5, -0.01))
    for a in ax:
        a.set_xticks([])
    save(fig, "fig3_alzette")


def fig4():
    p = os.path.join(WORK, "alzette_gauges.csv")
    if not os.path.exists(p):
        return
    g = pd.read_csv(p); g = g[g.output >= 1]
    fig, axs = plt.subplots(1, 5, figsize=(7.4, 2.2), constrained_layout=True, sharex=True)
    for a, name in zip(axs, GAUGES):
        for lab in STYLE:
            r = g[g.run == lab]
            if len(r):
                a.plot(r.output, r[name + "_Q"], **(STYLE[lab] if name == "Walferdange" else nolabel(STYLE[lab])))
        a.set_title(name, fontsize=8); a.set_xticks([0, 6, 12, 18, 24])
        a.set_xlabel("hours", fontsize=7); a.tick_params(labelsize=7)
    axs[0].set_ylabel("discharge (m$^3$ s$^{-1}$)")
    fig.legend(*axs[0].get_legend_handles_labels(), loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.12))
    save(fig, "fig4_gauges")


def fig5():
    cells = np.array([2e4, 8e4, 5e5, 2e6, 8e6]); res = ["10m", "5m", "2m", "1m", "0.5m"]

    def speed(lab):
        out = []
        for r, n in zip(res, cells):
            v = []                      # run times of the two repeats
            for rep in (1, 2):
                d = os.path.join(WORK, f"speed_{r}_{lab}_{rep}")
                try:
                    v.append((steps(d), loop_time(d)))
                except (OSError, IndexError, ValueError):
                    pass
            out.append(n * v[0][0] / np.median([t for _, t in v]) / 1e6 if v else np.nan)
        return np.array(out)

    s83, s82, sm = speed("cuda83"), speed("cuda82"), speed("metal")
    if np.isnan(s83).all() and np.isnan(sm).all():
        return
    fig, ax = plt.subplots(figsize=(4.4, 3.0), constrained_layout=True)
    ax.plot(cells, s83, "o-", color=CUDA, lw=1.6, ms=4, label="CUDA 8.3")
    ok = ~np.isnan(s82)
    if ok.any():
        ax.plot(cells[ok], s82[ok], "s--", color=CUDA, lw=1.0, ms=3.5, mfc="white", label="CUDA 8.2")
        if not ok.all():       # the first grid on which 8.2 stops
            i = int(np.argmin(ok))
            ax.plot(cells[i], s82[ok][-1], "x", color=CUDA, ms=7, mew=1.5)
            ax.annotate("8.2 stops\n(illegal address)", (cells[i], s82[ok][-1]), xytext=(6, -16),
                        textcoords="offset points", fontsize=7, color=CUDA)
    ax.plot(cells, sm, "o-", color=METAL, lw=1.6, ms=4, label="Metal")
    ax.set(xscale="log", yscale="log", xlabel="cells", ylabel="million cell updates per second")
    ax.legend(frameon=False, loc="upper left")
    save(fig, "fig5_speed")


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5):
        f()
