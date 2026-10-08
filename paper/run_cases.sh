#!/bin/sh
# Run the cases of the paper with one build.
#
#   paper/run_cases.sh BIN LABEL [CASES]
#
#   BIN    lisflood executable
#   LABEL  name of the platform, one of:
#            cpu     CPU FV1 solver           (no GPU)
#            cuda83  CUDA FV1, LISFLOOD-FP 8.3 (NVIDIA GPU)
#            cuda82  CUDA FV1, LISFLOOD-FP 8.2 (NVIDIA GPU; for the speed comparison)
#            metal   Metal FV1, LISFLOOD-FP 8.3 (Apple GPU)
#   CASES  comma list (default: dambreak,convergence,ea4,alzette,speed):
#            dambreak  Ritter dam break, 5 m, 60 s                  (seconds)
#            convergence  the dam break at 10, 5, 2.5, 1.25, 0.625 m  (seconds)
#            ea4       EA test 4 at 5 m and 1 m, 5 h                 (minutes on a GPU)
#            alzette   Alzette 10 m, 24 h; needs cases/alzette/ data (10 min to hours)
#            speed     EA test 4 at 10, 5, 2, 1 and 0.5 m, no output, twice each
#
# Results go to paper/work/<case>_<LABEL>/ (WORK= to change). Run it once per
# platform and collect the work/ folders in one place to compare them.
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1; LABEL=$2; CASES=${3:-dambreak,convergence,ea4,alzette,speed}
[ -x "$BIN" ] && [ -n "$LABEL" ] || { awk 'NR > 1 && /^#/ { print } /^HERE=/ { exit }' "$0"; exit 2; }
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")
WORK=${WORK:-$HERE/work}; mkdir -p "$WORK"
case $LABEL in
  cpu) GPU= ;;
  metal) GPU=-metal ;;
  *) GPU=-cuda ;;
esac
want() { echo ",$CASES," | grep -q ",$1,"; }

run() {  # name case-dir par
  O=$WORK/$1; rm -rf "$O"; mkdir -p "$O"
  case $2 in /*) C=$2 ;; *) C=$HERE/cases/$2 ;; esac
  (cd "$C" && "$BIN" $GPU -v -dirroot "$O" "$3" > "$O/run.log" 2>&1); rc=$?
  loop=$(awk '/^loop time/{print $3}' "$O/run.log" | tail -1)
  steps=$(tail -1 "$O"/*.mass 2>/dev/null | awk '{print $4}')
  echo "RESULT $1 rc=$rc loop_s=${loop:--} steps=${steps:--}"
}

want dambreak && run dambreak_$LABEL dambreak dambreak.par
if want convergence; then
  # the dam break of cases/dambreak (2000 m x 4 cells, dam at 1000 m) on finer grids,
  # depth written to 6 decimals; the time step is set by the Courant number alone
  for dx in 10 5 2.5 1.25 0.625; do
    C=$WORK/convergence_inputs/dx$dx; mkdir -p "$C"
    for f in dem start; do
      awk -v dx=$dx -v f=$f 'BEGIN { n = 2000 / dx
        printf "ncols %d\nnrows 4\nxllcorner 0.0\nyllcorner 0.0\ncellsize %s\nNODATA_value -9999\n", n, dx
        for (j = 0; j < 4; j++) { for (i = 0; i < n; i++) printf "%s ", (f == "start" && (i + 0.5) * dx < 1000) ? "1.0" : "0.0"; print "" } }' > "$C/$f.asc"
    done
    printf 'resroot r\nDEMfile dem.asc\nstartfile start.asc\nfv1\nfpfric 0.0\nsim_time 60\ninitial_tstep 10\nmassint 10\nsaveint 60\noutput_precision 6\n' > "$C/db.par"
    run convergence_${dx}_$LABEL "$C" db.par
  done
fi
if want ea4 || want speed; then
  [ -f "$HERE/cases/ea_test4/ea4-5m.dem" ] || "$HERE/cases/ea_test4/get_data.sh"
fi
if want ea4; then
  run ea4_5m_$LABEL ea_test4 ea4_5m.par
  run ea4_1m_$LABEL ea_test4 ea4_1m.par
fi
if want alzette; then
  A=$HERE/cases/alzette
  [ -f "$A/Alzette_sub_basin_10m_bridge_burn.asc" ] || [ ! -f "$A/alzette_inputs.tar.gz" ] ||
    (cd "$A" && tar -xzf alzette_inputs.tar.gz && echo "alzette: inputs unpacked")
  if [ -f "$HERE/cases/alzette/Alzette_sub_basin_10m_bridge_burn.asc" ]; then
    run alzette_$LABEL alzette alzette.par
  else
    echo "alzette: input files not found in cases/alzette (see cases/alzette/README.md); skipped"
  fi
fi
if want speed; then
  for r in 10m 5m 2m 1m 0.5m; do
    for rep in 1 2; do run speed_${r}_${LABEL}_$rep ea_test4 speed_$r.par; done
  done
fi
