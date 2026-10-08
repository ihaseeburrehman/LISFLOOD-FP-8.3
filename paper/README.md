# Reproduce the paper

The three cases of the paper, run with the FV1 solver on the CPU, an NVIDIA GPU
(CUDA) and an Apple GPU (Metal), and the scripts that compute every number of
its Sects. 4 and 5 and draw its Figs. 1 to 5.

```
run_cases.sh          runs the cases with one build
cases/dambreak/       Ritter dam break, 5 m, 60 s (inputs included); run_cases.sh also
                      runs it at 10 to 0.625 m for the convergence test
cases/ea_test4/       EA test 4 at 5 m and 1 m, and the speed runs at 10 m to 0.5 m
                      (get_data.sh downloads the DEMs and the published reference)
cases/alzette/        July 2021 Alzette flood, 10 m, 24 h (inputs included, packed)
analysis/results.py   prints every number of the paper
analysis/figures.py   draws Figs. 1 to 5
analysis/common.py    shared helpers
```

## 1. Run the cases on each platform

Build LISFLOOD-FP 8.3 as in the top-level README, then on each machine
(`LABEL` names the platform: use these four, so that runs from different
machines can be put side by side):

```bash
paper/run_cases.sh build/lisflood cpu        # CPU FV1, any machine
```

```bash
paper/run_cases.sh build/lisflood cuda83     # NVIDIA GPU, LISFLOOD-FP 8.3
```

```bash
paper/run_cases.sh build82/lisflood cuda82 speed   # NVIDIA GPU, released 8.2 (speed comparison only)
```

```bash
paper/run_cases.sh build/lisflood metal      # Apple GPU
```

To build the released 8.2 for comparison with CUDA 12, rename one unused
function first (8.2 does not compile otherwise; nothing else changes), from the
top folder of the unpacked 8.2:

```bash
sed -i.bak 's/cudaError_t sync()/cudaError_t sync_device()/' cuda/adaptive/utilities/cuda_utils.cu cuda/adaptive/utilities/cuda_utils.cuh
```

A third argument limits the cases: `dambreak`, `convergence`, `ea4`,
`alzette`, `speed` (comma separated). Results go to `paper/work/<case>_<label>/`. The dam break
takes seconds; EA test 4 minutes on a GPU (an hour on a CPU at 1 m); the
Alzette day about 10 min on an NVIDIA A100, 33 min on an Apple M1 Max and 4 h
on 128 CPU cores; the speed runs about 10 min on an A100.

## 2. Collect and analyse

Copy the `work/` folders of all machines into one `paper/work/`, then (Python 3
with numpy, pandas and matplotlib):

```bash
python3 paper/analysis/results.py
```

```bash
python3 paper/analysis/figures.py
```

`results.py` prints the numbers of Table 4 and Sects. 4 and 5; `figures.py`
writes `fig1_dambreak` to `fig5_speed` (pdf and png) into `paper/work/`. A run
that is missing is skipped. Run times depend on the machine; the other numbers
should match the paper, up to differences of single-precision rounding on a
different GPU model.
