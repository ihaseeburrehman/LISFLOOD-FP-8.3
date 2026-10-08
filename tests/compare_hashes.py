#!/usr/bin/env python3
"""Compare runs of a test suite by the fingerprints from hash_outputs.sh.

    compare_hashes.py hashes.txt REF LABEL [LABEL ...] [--cases A,B] [--md FILE]

For every <case>_<solver> that has a REF run, each LABEL run is reported as

  bit-identical      every output file has the same fingerprint
  same except .mass  only the .mass log differs (its flux columns are sums of
                     floating-point atomics, whose last digit varies)
  differs            other files differ; they are listed, and compare_suite.py
                     gives the size of the differences for those runs
  different files    the two runs did not write the same set of files

Only fingerprints are read, so this takes seconds.
"""
import sys
from collections import defaultdict


def main():
    args = sys.argv[1:]
    md, only = None, None
    for opt in ("--md", "--cases"):
        if opt in args:
            i = args.index(opt)
            val = args[i + 1]
            del args[i:i + 2]
            if opt == "--md":
                md = val
            else:
                only = set(val.split(","))
    hashes, ref, labels = args[0], args[1], args[2:]

    runs = defaultdict(dict)          # (case_solver) -> label -> {file: md5}
    with open(hashes) as f:
        for line in f:
            md5, path = line.split(None, 1)
            run, _, name = path.strip().partition("/")
            key, _, lab = run.rpartition("_")
            runs[key].setdefault(lab, {})[name] = md5

    rows = ["| case | solver | run | result | files that differ |", "|---|---|---|---|---|"]
    for key in sorted(runs):
        case, _, solver = key.partition("_")
        if only and case not in only or ref not in runs[key]:
            continue
        r = runs[key][ref]
        for lab in labels:
            o = runs[key].get(lab)
            if o is None:
                continue
            if set(o) != set(r):
                res, detail = "different files", ""
            else:
                diff = sorted(f for f in r if r[f] != o[f])
                if not diff:
                    res, detail = "bit-identical", ""
                elif all(f.endswith(".mass") for f in diff):
                    res, detail = "same except .mass", ""
                else:
                    res = "differs"
                    detail = ", ".join(diff[:6]) + (" ..." if len(diff) > 6 else "")
            rows.append("| %s | %s | %s vs %s | %s | %s |" % (case, solver, lab, ref, res, detail))
    text = "\n".join(rows)
    print(text)
    if md:
        with open(md, "w") as f:
            f.write(text + "\n")


if __name__ == "__main__":
    main()
