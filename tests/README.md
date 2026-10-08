# Tests

Run LISFLOOD-FP's own test cases (from `testing/` of the 8.2 source, in
`data/`) with any build.

```
quick_test.sh      a few minutes: does my build start, run and finish as expected?
run_suite.sh       every case with every solver that applies (list: cases.txt)
hash_outputs.sh    fingerprints every output file of a suite, once, in parallel
compare_hashes.py  compares two builds file by file from the fingerprints
compare_suite.py   size of the differences (mm) where files differ
refine_dem.py      makes the 5.2 million-cell case T07x3 (used by run_suite.sh)
cases.txt          the cases, their solvers and how their inputs are prepared
data/              the test inputs
```

## Quick test

```bash
tests/quick_test.sh build/lisflood              # NVIDIA GPU: FV1, FV2, DG2, ACC
```

```bash
GPU=metal tests/quick_test.sh build/lisflood    # Apple GPU: FV1
```

It runs four small cases (fv1, T01, T03, T031) with every solver of the back
end and checks that each run finishes and takes the number of time steps of the
verified runs (NVIDIA A100, Apple M1 Max), within 1 %. It is a check that the
build works, not a comparison of results: for that, use the suite below.

## Compare two builds (e.g. 8.2 and 8.3) on every case

To build the released 8.2 for comparison with CUDA 12, rename one unused
function first (8.2 does not compile otherwise; nothing else changes), from the
top folder of the unpacked 8.2:

```bash
sed -i.bak 's/cudaError_t sync()/cudaError_t sync_device()/' cuda/adaptive/utilities/cuda_utils.cu cuda/adaptive/utilities/cuda_utils.cuh
```

```bash
tests/run_suite.sh build82/lisflood S1 && tests/run_suite.sh build/lisflood N1 && tests/run_suite.sh build/lisflood N2
```

```bash
tests/hash_outputs.sh work/out > hashes.txt && python3 tests/compare_hashes.py hashes.txt S1 N1 N2
```

Each run is reported as `bit-identical`, `same except .mass` (only the
mass-balance log differs: its flux columns are sums that many GPU threads add
to at once, so their last digit varies between any two runs), `differs` or
`different files` (e.g. a run that stopped). N1 against N2 checks that 8.3
gives identical results in repeated runs. Cases T05, T06, T07 and T07x3 are
prepared by `generateDG2DEM`, whose output can differ slightly between
operating systems; to compare two machines, prepare once
(`tests/run_suite.sh build/lisflood prep none`) and copy `work/data`.
