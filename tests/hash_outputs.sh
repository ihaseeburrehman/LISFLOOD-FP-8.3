#!/bin/sh
# Fingerprint every output file of a test suite once, in parallel.
#
#   hash_outputs.sh WORK/out [JOBS] > hashes.txt
#
# Writes one line per file: "<md5>  <run-directory>/<file>". run.log and run.par
# are left out. Each file is read once, however many runs it is compared with,
# so this is much faster than comparing runs pair by pair; compare_hashes.py
# then compares the fingerprints. JOBS parallel readers (default 32).
OUT=$1; JOBS=${2:-32}
[ -d "$OUT" ] || { sed -n '2,10p' "$0"; exit 2; }
cd "$OUT" && find . -mindepth 2 -maxdepth 2 -type f ! -name run.log ! -name run.par -print0 |
  xargs -0 -n 16 -P "$JOBS" md5sum | sed 's#  \./#  #' | sort -k2
