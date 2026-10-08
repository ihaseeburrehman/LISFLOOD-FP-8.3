#!/usr/bin/env python3
"""Print every number the paper reports for its three cases and run times.

    python3 paper/analysis/results.py [WORK]

WORK (default paper/work) holds the run folders made by run_cases.sh on each
platform: <case>_<label>, with labels cpu, cuda82, cuda83 and metal. Missing
runs are skipped. Writes WORK/alzette_gauges.csv (used by figures.py).
"""
import os
import sys

import numpy as np
import pandas as pd

from common import CASES, GAUGES, grid, gauges, header, loop_time, outputs, ritter, steps

WORK = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "work"))
LABELS = ["cpu", "cuda82", "cuda83", "metal"]


def run(case, label):
    d = os.path.join(WORK, f"{case}_{label}")
    return d if os.path.isdir(d) and outputs(d, "wd") else None


def stats(a, b, wet):
    d = (a - b)[wet]
    return 100 * (d == 0).mean(), 1000 * np.sqrt((d ** 2).mean()), 1000 * np.abs(d).max()


def dambreak():
    print("\n== Dam break (5 m, 400 x 4 cells): RMSE to the exact solution, mm")
    x = (np.arange(400) + 0.5) * 5.0
    depth = {}
    for lab in LABELS:
        d = run("dambreak", lab)
        if not d:
            continue
        g = [grid(outputs(d, "wd")[k]) for k in (1, 2, 3)]
        depth[lab] = g
        rmse = [1000 * np.sqrt(((g[i][1] - ritter(x, t)) ** 2).mean()) for i, t in enumerate((20, 40, 60))]
        print(f"  {lab:7s} steps {steps(d):4d}   20 s {rmse[0]:.2f}   40 s {rmse[1]:.2f}   60 s {rmse[2]:.2f}")
    labs = list(depth)
    for a in labs[1:]:
        same = all(np.array_equal(depth[labs[0]][i], depth[a][i]) for i in range(3))
        print(f"  {a} and {labs[0]}: same depth in every cell at 20, 40 and 60 s: {same}")


def convergence():
    print("\n== Dam break convergence at 60 s: error against the exact solution, mm, and observed order")
    first = {}
    for lab in LABELS:
        prev = None
        for dx in ("10", "5", "2.5", "1.25", "0.625"):
            d = run(f"convergence_{dx}", lab)
            if not d:
                continue
            h = grid(outputs(d, "wd")[1]).astype(np.float64)
            x = (np.arange(h.shape[1]) + 0.5) * float(dx)
            e = h[1] - ritter(x, 60.0)
            l1, l2 = 1000 * np.abs(e).mean(), 1000 * np.sqrt((e ** 2).mean())
            o = "" if prev is None else f"   order L1 {np.log2(prev[0] / l1):.2f}  RMSE {np.log2(prev[1] / l2):.2f}"
            print(f"  {lab:7s} dx {dx:>5s} m  steps {steps(d):5d}  L1 {l1:.3f}  RMSE {l2:.3f}{o}")
            prev = (l1, l2)
            if lab != "cpu" and run(f"convergence_{dx}", "cpu"):
                c = grid(outputs(run(f"convergence_{dx}", "cpu"), "wd")[1])
                first.setdefault(lab, []).append(1000 * float(np.abs(grid(outputs(d, "wd")[1]) - c).max()))
    for lab, v in first.items():
        print(f"  {lab} against cpu: largest depth difference {max(v):.4f} mm")


def ea4():
    for res in ("5m", "1m"):
        ref = os.path.join(CASES, "ea_test4", f"reference-{res}")
        if not os.path.isdir(ref):
            print(f"\n== EA test 4 {res}: published reference missing (cases/ea_test4/get_data.sh)")
            continue
        print(f"\n== EA test 4, {res}: against the published run (CPU FV1 of 8.0, double precision)")
        print("   published steps:", steps(ref))
        for lab in LABELS:
            d = run(f"ea4_{res}", lab)
            if not d:
                continue
            rm, mx = [], []
            for k in range(1, 6):
                a, b = grid(outputs(d, "wd")[k]), grid(outputs(ref, "wd")[k])
                _, r, m = stats(a, b, (a > 0.001) | (b > 0.001))
                rm.append(r); mx.append(m)
            print(f"  {lab:7s} steps {steps(d):7d}   RMSE (worst hour) {max(rm):.2f} mm   largest difference {max(mx):.0f} mm")


def alzette():
    cpu = run("alzette", "cpu")
    if not cpu:
        print("\n== Alzette: no CPU run, skipped")
        return
    print("\n== Alzette 10 m, 24 h: against the CPU, 24 hourly outputs, wet cells (> 1 cm)")
    dem = os.path.join(CASES, "alzette", "Alzette_sub_basin_10m_bridge_burn.asc")
    if not os.path.isfile(dem):
        print("   input files not found in cases/alzette (see cases/alzette/README.md); skipped")
        return
    basin = grid(dem) > -9000
    print(f"   cells {basin.size:,}, in the basin {int(basin.sum()):,}; CPU steps {steps(cpu)}")
    rows = []
    for lab in LABELS:
        d = run("alzette", lab)
        if not d:
            continue
        rows.append(gauges(d).assign(run=lab))
        if lab == "cpu":
            continue
        n = same = 0; sq = 0.0; mx = 0.0; hourly = []
        for k in range(1, 25):
            a, b = grid(outputs(d, "wd")[k]), grid(outputs(cpu, "wd")[k])
            wet = basin & ((a > 0.01) | (b > 0.01))
            dd = (a - b)[wet]
            n += dd.size; same += int((dd == 0).sum()); sq += float((dd.astype(np.float64) ** 2).sum())
            mx = max(mx, float(np.abs(dd).max())); hourly.append(100 * (dd == 0).mean())
            if k == 7:
                diff = np.abs(a - b)[wet] * 1000
                peak = (int(wet.sum()), int((diff > 0).sum()), int(((diff > 0) & (diff <= 1.5)).sum()),
                        int(((diff > 1.5) & (diff <= 10.5)).sum()), int((diff > 10.5).sum()))
        print(f"  {lab:7s} steps {steps(d)}   identical {100 * same / n:.1f} % "
              f"(hours {min(hourly):.1f}-{max(hourly):.1f} %)   RMSE {1000 * np.sqrt(sq / n):.2f} mm   "
              f"largest {1000 * mx:.0f} mm")
        print(f"          at 15 July 01:00: {peak[1]:,} of {peak[0]:,} wet cells differ: "
              f"{peak[2]:,} by 1 mm, {peak[3]} by 2-10 mm, {peak[4]} by more")
    g = pd.concat(rows)
    g.to_csv(os.path.join(WORK, "alzette_gauges.csv"), index=False)
    c = g[(g.run == "cpu") & (g.output >= 1)].set_index("output")
    print("   gauges, against the CPU: peak discharge difference, largest |dQ| (m3/s), largest |dh| (mm)")
    for lab in [l for l in LABELS if l != "cpu" and (g.run == l).any()]:
        r = g[(g.run == lab) & (g.output >= 1)].set_index("output")
        for n in GAUGES:
            pc, pr = c[n + "_Q"].max(), r[n + "_Q"].max()
            print(f"  {lab:7s} {n:12s} peak {pc:7.2f} vs {pr:7.2f} ({100 * (pr - pc) / pc:+.3f} %), "
                  f"same hour {c[n + '_Q'].idxmax() == r[n + '_Q'].idxmax()}, "
                  f"|dQ| {(r[n + '_Q'] - c[n + '_Q']).abs().max():.3f}, "
                  f"|dh| {1000 * (r[n + '_depth'] - c[n + '_depth']).abs().max():.0f}")


def volume(path):
    """Water volume stored in a depth grid, m3."""
    a = grid(path).astype(np.float64)
    return float(a[a > 0].sum()) * header(path)["cellsize"] ** 2


def ea4_inflow(t):
    """Water added through the breach of EA test 4 by time t (ea4.bdy: 20 m wide,
    0 until 300 s, rising linearly to 1 m2/s at 3600 s, then constant)."""
    return 0.0 if t < 300 else (20 * (t - 300) ** 2 / 6600 if t < 3600 else 33000 + 20 * (t - 3600))


def mass():
    print("\n== Volume in the depth grids against the water present or added, %")
    for lab in LABELS:
        d = run("dambreak", lab)
        if d:
            print(f"  dam break  {lab:7s} " + "  ".join(
                f"{20 * k} s {100 * (volume(outputs(d, 'wd')[k]) - 2e4) / 2e4:+.4f}" for k in (1, 2, 3)))
    for res in ("5m", "1m"):
        for lab in LABELS:
            d = run(f"ea4_{res}", lab)
            if d:
                print(f"  EA4 {res:3s}    {lab:7s} " + "  ".join(
                    f"{k} h {100 * (volume(outputs(d, 'wd')[k]) / ea4_inflow(3600 * k) - 1):+.4f}" for k in (1, 2)))


def run_times():
    print("\n== Run times (time loop, s)")
    for lab in LABELS:
        d = os.path.join(WORK, f"alzette_{lab}")
        if os.path.isdir(d):
            t = loop_time(d)
            print(f"  Alzette {lab:7s} {t:8.0f} s = {t / 60:6.1f} min")
    print("  EA test 4 without output: million cell updates per second (average run time of two runs)")
    cells = {"10m": 2e4, "5m": 8e4, "2m": 5e5, "1m": 2e6, "0.5m": 8e6}
    for lab in LABELS:
        line, spread = [], []
        for r, n in cells.items():
            v = []                      # run times of the two repeats
            for rep in (1, 2):
                d = os.path.join(WORK, f"speed_{r}_{lab}_{rep}")
                try:
                    v.append((steps(d), loop_time(d)))
                except (OSError, IndexError, ValueError):
                    pass
            line.append(f"{r} {n * v[0][0] / np.mean([t for _, t in v]) / 1e6:6.0f}" if v else f"{r}   stops")
            if len(v) == 2 and n >= 5e5:
                spread.append(f"{r} {abs(v[0][1] - v[1][1]):.0f} s ({100 * abs(v[0][1] - v[1][1]) / min(t for _, t in v):.1f} %)")
        if any("stops" not in s for s in line):
            print(f"  {lab:7s} " + "  ".join(line))
            print(f"          difference between the two runs: " + "  ".join(spread))


if __name__ == "__main__":
    print("run folders in", WORK)
    dambreak(); convergence(); ea4(); mass(); alzette(); run_times()
