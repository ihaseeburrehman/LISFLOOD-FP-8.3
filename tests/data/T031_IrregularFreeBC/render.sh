#!/bin/bash
for f in results/*.elev; do
	gdal_translate -of png $f $f.png -outsize 1600% 1600% -scale 0 0.2 0 255
done
