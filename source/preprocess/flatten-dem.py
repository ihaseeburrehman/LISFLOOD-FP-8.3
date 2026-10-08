#!/usr/bin/env python3
from osgeo import gdal
import math
import sys
import numpy as np

def main():
    if len(sys.argv) <= 3:
        print('Usage: speed.py <dem_in> <dem_out> <flattened_value>')
        sys.exit(1)

    gdal.SetConfigOption('AAIGRID_DATATYPE', 'Float64')
    in_file = gdal.Open(sys.argv[1], gdal.GA_ReadOnly)
    in_band = in_file.GetRasterBand(1)
    in_field = in_band.ReadAsArray()

    transform = in_file.GetGeoTransform()
    nodata = in_band.GetNoDataValue()

    new_value = sys.argv[3]

    with open(sys.argv[2], 'w') as out:
        print("ncols", in_field.shape[1], file=out)
        print("nrows", in_field.shape[0], file=out)
        print("xllcorner", transform[0], file=out)
        print("yllcorner", transform[3]-transform[1]*in_field.shape[0], file=out)
        print("cellsize", transform[1], file=out)
        print("NODATA_value", nodata, file=out)

        for in_row in in_field:
            for old_value in in_row:
                if (old_value == nodata):
                    print(nodata, end=' ', file=out)
                else:
                    print(new_value, end=' ', file=out)
            print(file=out)

    in_file = None

if __name__ == '__main__':
    main()
