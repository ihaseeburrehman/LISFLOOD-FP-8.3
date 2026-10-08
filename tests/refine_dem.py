#!/usr/bin/env python3
"""Split every cell of an ASCII raster into N x N cells of the same value.

    refine_dem.py RASTER N

The raster is overwritten. The terrain is unchanged, only sampled finer, so the
boundary and gauge files of the case still apply. Used to turn LISFLOOD-FP's
Carlisle test (951 x 611 cells) into a grid large enough to exercise the
large-grid code paths.
"""
import sys


def main():
    path, n = sys.argv[1], int(sys.argv[2])
    with open(path) as f:
        head = [next(f) for _ in range(6)]
        rows = [line.split() for line in f if line.strip()]
    key = [h.split()[0] for h in head]
    val = [h.split()[1] for h in head]
    ncols, nrows = int(val[0]), int(val[1])
    assert len(rows) == nrows and all(len(r) == ncols for r in rows), "unexpected raster shape"
    val[0], val[1] = str(ncols * n), str(nrows * n)
    val[4] = repr(float(val[4]) / n)
    with open(path, "w") as f:
        for k, v in zip(key, val):
            f.write("%-14s %s\n" % (k, v))
        for r in rows:
            line = " ".join(c for c in r for _ in range(n)) + "\n"
            for _ in range(n):
                f.write(line)


if __name__ == "__main__":
    main()
