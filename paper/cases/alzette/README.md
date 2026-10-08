# Alzette, July 2021, 10 m, 24 h

The upper Alzette basin (Luxembourg) on a 10 m grid of 2440 x 4362 cells, from
14 July 2021 18:00 to 15 July 18:00 UTC, forced by rainfall from the Weather
Research and Forecasting model. `alzette.par` and `Border.bci` are the run
settings of the paper.

The input files are in `alzette_inputs.tar.gz` (71 MB; 260 MB unpacked), which
`run_cases.sh` unpacks into this folder the first time the case runs:

```
Alzette_sub_basin_10m_bridge_burn.asc   terrain, bridges removed
manning.n.ascii                         Manning's n by land use
rain_wrf1p3_24h_from_14jul18.nc         rainfall, NetCDF
start_14jul_18hr.wd                     initial water depth
```

They are the inputs of the runs in the paper. How they were made is described
in the companion studies cited in the paper.
