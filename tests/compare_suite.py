#!/usr/bin/env python3
"""Compare the runs made by run_suite.sh.

    compare_suite.py WORK/out REF LABEL [LABEL ...]  [--md FILE] [--cases A,B]

For every <case>_<solver> that has a REF run, each LABEL run is compared with it:

  bit-identical      every output file is byte-for-byte the same
  same except .mass  only the .mass log differs (its flux and error columns are
                     sums of floating-point atomics, whose order varies)
  differs            then: share of identical cells, largest difference and RMSE
                     over wet cells (> 1 mm) for the last depth grid and .max

run.log and run.par are not compared. Results print as a Markdown table, and are
also written to FILE with --md. --cases limits the comparison to those cases.
"""
import os
import sys

import numpy as np

SKIP = {"run.log", "run.par"}


def read_grid(path):
    with open(path) as f:
        head = [next(f) for _ in range(6)]
    nd = float(head[5].split()[1]) if head[5].lower().startswith("nodata") else -9999.0
    return np.loadtxt(path, skiprows=6, dtype=np.float64), nd


def grid_stats(a_path, b_path):
    a, na = read_grid(a_path)
    b, nb = read_grid(b_path)
    if a.shape != b.shape:
        return "shape differs"
    m = (a != na) & (b != nb)
    d = np.abs(a - b)[m]
    wet = m & ((a > 0.001) | (b > 0.001))
    dw = (a - b)[wet]
    rmse = np.sqrt((dw ** 2).mean()) if dw.size else 0.0
    return "%.2f%% same, max %.3g mm, RMSE %.3g mm" % (100 * (d == 0).mean(), 1000 * d.max(), 1000 * rmse)


def steps(d):
    for f in os.listdir(d):
        if f.endswith(".mass"):
            with open(os.path.join(d, f)) as fh:
                lines = fh.read().strip().splitlines()
            return lines[-1].split()[3] if len(lines) > 1 else "-"
    return "-"


def compare(rdir, odir):
    rfiles = sorted(f for f in os.listdir(rdir) if f not in SKIP)
    ofiles = sorted(f for f in os.listdir(odir) if f not in SKIP)
    if not rfiles:
        return "reference has no output", ""
    if rfiles != ofiles:
        return "different file sets", ""
    try:
        diff = [f for f in rfiles
                if open(os.path.join(rdir, f), "rb").read() != open(os.path.join(odir, f), "rb").read()]
    except OSError as e:  # e.g. a run still being written, or a file-system error
        return "not readable", str(e)
    if not diff:
        return "bit-identical", ""
    if all(f.endswith(".mass") for f in diff):
        return "same except .mass", ""
    wd = [f for f in rfiles if f.endswith(".wd")]
    detail = []
    for f in ([wd[-1]] if wd else []) + [f for f in rfiles if f.endswith(".max")]:
        detail.append("%s: %s" % (f, grid_stats(os.path.join(rdir, f), os.path.join(odir, f))))
    return "differs (%d files)" % len(diff), "; ".join(detail)


def main():
    args = sys.argv[1:]
    md, only = None, None
    if "--md" in args:
        i = args.index("--md")
        md = args[i + 1]
        del args[i:i + 2]
    if "--cases" in args:
        i = args.index("--cases")
        only = set(args[i + 1].split(","))
        del args[i:i + 2]
    out, ref, labels = args[0], args[1], args[2:]
    runs = {}
    for d in sorted(os.listdir(out)):
        key, _, lab = d.rpartition("_")
        runs.setdefault(key, {})[lab] = os.path.join(out, d)
    rows = ["| case | solver | %s steps | run | steps | result | detail |" % ref,
            "|---|---|---|---|---|---|---|"]
    for key in sorted(runs):
        if ref not in runs[key]:
            continue
        case, _, solver = key.partition("_")
        if only and case not in only:
            continue
        r = runs[key][ref]
        for lab in labels:
            if lab not in runs[key]:
                continue
            res, detail = compare(r, runs[key][lab])
            rows.append("| %s | %s | %s | %s | %s | %s | %s |" % (
                case, solver, steps(r), lab, steps(runs[key][lab]), res, detail))
    text = "\n".join(rows)
    print(text)
    if md:
        with open(md, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
