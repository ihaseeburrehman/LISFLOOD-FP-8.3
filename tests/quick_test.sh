#!/bin/sh
# Quick check of a fresh LISFLOOD-FP 8.3 build: a few minutes on one GPU.
#
#   tests/quick_test.sh build/lisflood          # NVIDIA GPU (CUDA): FV1, FV2, DG2, ACC
#   GPU=metal tests/quick_test.sh build/lisflood   # Apple GPU (Metal): FV1
#
# Prints the version banner, runs four small test cases that come with
# LISFLOOD-FP (fv1, T01, T03, T031) with every solver of the back end, and
# checks that each run used the GPU and finished, and that it took the number of
# time steps of the verified runs (NVIDIA A100, Apple M1 Max) within 1 %; a
# different GPU model may differ by a few steps in single precision.
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1; GPU=${GPU:-cuda}
[ -x "$BIN" ] || { sed -n '2,12p' "$0"; exit 2; }
"$BIN" -version | grep 'LISFLOOD-FP'

if [ "$GPU" = metal ]; then
  SOLVERS=fv1; GPU_USED="Metal: Using GPU"
  EXPECT="fv1_fv1=10000 T01_fv1=3517 T03_fv1=532 T031_fv1=206"
else
  SOLVERS=all; GPU_USED="Using CUDA device"
  EXPECT="fv1_fv1=10000 fv1_fv2=10000 fv1_acc=10000
          T01_fv1=3516 T01_fv2=5327 T01_dg2=5327 T01_acc=2471
          T03_fv1=531 T03_fv2=798 T03_dg2=912 T03_acc=4929
          T031_fv1=211 T031_fv2=219 T031_dg2=219 T031_acc=200"
fi
WORK=${WORK:-$PWD/quick_test_work}; export WORK GPU
"$HERE/run_suite.sh" "$BIN" Q "$SOLVERS" fv1,T01,T03,T031 > "$WORK.log" 2>&1

fail=0
for e in $EXPECT; do
  key=${e%=*}; want=${e#*=}
  got=$(awk -v k="$key" '{split($2,a," ")} $2"_"$3==k {print $7}' "$WORK.log" | sed 's/steps=//')
  rc=$(awk -v k="$key" '$2"_"$3==k {print $5}' "$WORK.log")
  ok=$(awk -v g="$got" -v w="$want" 'BEGIN{print (g!="" && g!="-" && (g-w)^2 <= (0.01*w)^2) ? "ok" : "FAIL"}')
  [ "$rc" = "rc=0" ] || ok=FAIL
  grep -q "$GPU_USED" "$WORK/out/${key}_Q/run.log" 2>/dev/null || ok="FAIL (GPU not used)"
  [ "$ok" = ok ] || fail=1
  printf "  %-8s %-4s steps %8s  (verified %s)  %s\n" "${key%_*}" "${key#*_}" "${got:--}" "$want" "$ok"
done
[ $fail = 0 ] && echo "QUICK TEST PASSED" || echo "QUICK TEST FAILED: see $WORK.log and $WORK/out/*/run.log"
exit $fail
