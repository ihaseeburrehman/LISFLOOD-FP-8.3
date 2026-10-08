#!/bin/sh
# Download EA test 4 (Neelz and Pender 2013; Shaw et al. 2021, Sect. 3.1) from
# Zenodo (doi 10.5281/zenodo.4661888, CC-BY-4.0, about 17 MB): the DEMs at
# 10, 5, 2, 1 and 0.5 m, and the published LISFLOOD-FP 8.0 output of the CPU FV1
# solver at 5 m and 1 m, which is the reference of the paper.
set -e
cd "$(dirname "$0")"
[ -f ea4.zip ] || curl -L -o ea4.zip https://zenodo.org/records/4661888/files/ea4.zip
unzip -q -o ea4.zip
for r in 10m 5m 2m 1m 0.5m; do gzip -dc "ea4/ea4-$r.dem.gz" > "ea4-$r.dem"; done
rm -rf reference-5m reference-1m
cp -r ea4/ea4-5m/fv1 reference-5m
cp -r ea4/ea4-1m/fv1 reference-1m
rm -rf ea4 ea4.zip
echo "EA test 4: DEMs and published reference (reference-5m, reference-1m) ready"
