#!/bin/sh
# Run LISFLOOD-FP's own test cases (see cases.txt) with one executable.
#
#   run_suite.sh BIN LABEL [SOLVERS] [CASES]
#
#   BIN      lisflood executable; generateDG2DEM/generateDG2start are taken from
#            the same directory
#   LABEL    names this run in the output, e.g. S1 (8.2, first run), N2 (8.3,
#            second run), M1 (Metal), C1 (CPU)
#   SOLVERS  comma list to restrict the solvers, e.g. fv1; "all" or empty: per case
#   CASES    comma list to restrict the cases, e.g. T01,T03; "all" or empty: every case
#
#   TESTING  folder holding the test inputs (default: data/ next to this script)
#   WORK     where prepared inputs and results go (default: ./work)
#   GPU      cuda (default), metal, or none for the CPU solvers
#   WRAP     command to run lisflood under, e.g. "compute-sanitizer --tool memcheck"
#   EXTRA    extra command-line options for lisflood, e.g. -metal (with GPU=none)
#
# Each run writes WORK/out/<case>_<solver>_<LABEL>/ and prints one RESULT line.
# SOLVERS=none prepares the inputs of the selected cases and runs nothing; do
# that once before starting several runs in parallel.
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$1; LABEL=$2; ONLY_S=$3; ONLY_C=$4
[ -n "$BIN" ] && [ -n "$LABEL" ] || { sed -n '2,22p' "$0"; exit 2; }
TESTING=${TESTING:-$HERE/data}
BIN=$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")
TOOLS=$(dirname "$BIN")
WORK=${WORK:-$PWD/work}; GPU=${GPU:-cuda}
mkdir -p "$WORK/data" "$WORK/out"

keep() { [ -z "$1" ] || [ "$1" = all ] || echo ",$1," | grep -q ",$2,"; }
par_value() { awk -v k="$2" 'tolower($1)==tolower(k){print $2}' "$1" | tr -d '\r'; }

prepare() {  # case dir par prep
  D=$WORK/data/$1
  [ -f "$D/.prepared" ] && return 0
  rm -rf "$D"; mkdir -p "$D"; cp -R "$TESTING/$2/." "$D/"
  for g in "$D"/*.gz; do [ -e "$g" ] && gunzip -f "$g"; done
  case $4 in *refine3*)
    raw="$D/$(par_value "$D/$3" DEMfile).raw"
    python3 "$HERE/refine_dem.py" "$raw" 3 || { echo "refine failed for $1" >&2; return 1; } ;;
  esac
  case $4 in *dem*)
    (cd "$D" && "$TOOLS/generateDG2DEM" "$3" > prepare_dem.log 2>&1) || { echo "generateDG2DEM failed for $1" >&2; return 1; } ;;
  esac
  case $4 in *start*)
    (cd "$D" && "$TOOLS/generateDG2start" "$3" > prepare_start.log 2>&1) || { echo "generateDG2start failed for $1" >&2; return 1; } ;;
  esac
  touch "$D/.prepared"
}

make_par() {  # src-par fixes solver -> stdout
  awk -v fixes="$2" -v solver="$3" -v gpu="$GPU" '
    BEGIN { n = split(fixes, f, ",")
            for (i = 1; i <= n; i++) {
              if (f[i] ~ /^drop:/) { sub(/^drop:/, "", f[i]); drop[tolower(f[i])] = 1 }
              else if (f[i] ~ /=/) { split(f[i], kv, "="); set[tolower(kv[1])] = kv[2] }
              else if (f[i] != "-") drop[tolower(f[i])] = 1 } }
    { sub(/\r$/, ""); k = tolower($1) }
    k == "dirroot" || k == "fv1" || k == "fv2" || k == "dg2" || k == "acceleration" ||
    k == "cuda" || k == "metal" || (k in drop) { next }
    (k in set) { print $1 "\t" set[k]; next }
    { print }
    END { print (solver == "acc" ? "acceleration" : solver)
          if (gpu != "none") print gpu }' "$1"
}

grep -v '^#' "$HERE/cases.txt" | while read -r name dir par solvers prep fixes; do
  [ -n "$name" ] || continue
  keep "$ONLY_C" "$name" || continue
  prepare "$name" "$dir" "$par" "$prep" || continue
  for s in $(echo "$solvers" | tr , ' '); do
    keep "$ONLY_S" "$s" || continue
    O=$WORK/out/${name}_${s}_$LABEL; rm -rf "$O"; mkdir -p "$O"
    make_par "$WORK/data/$name/$par" "${fixes:--}" "$s" > "$O/run.par"
    start=$(date +%s)
    (cd "$WORK/data/$name" && $WRAP "$BIN" $EXTRA -v -dirroot "$O" "$O/run.par" > "$O/run.log" 2>&1); rc=$?
    wall=$(( $(date +%s) - start ))
    steps=$(tail -1 "$O"/*.mass 2>/dev/null | awk '{print $4}')
    echo "RESULT $name $s $LABEL rc=$rc wall_s=$wall steps=${steps:--}"
  done
done
