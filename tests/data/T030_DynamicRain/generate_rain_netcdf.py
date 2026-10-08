#!/usr/bin/env python3
import netCDF4
import numpy as np

dt = 1/3600 # hours
dy = dx = 3 # metres
blx = 0
bly = 0
xsz = 2
ysz = 4
time_levels = 3

outfile = netCDF4.Dataset('dynamicrain.nc', 'w')
outfile.createDimension(dimname='time', size=time_levels)
xdim = outfile.createDimension('x', xsz)
ydim = outfile.createDimension('y', ysz)

times = outfile.createVariable('time', datatype='float32', dimensions=('time'))
xs = outfile.createVariable('x', datatype='float64', dimensions=('x'))
ys = outfile.createVariable('y', datatype='float64', dimensions=('y'))
rainfall = outfile.createVariable('rainfall_depth', dimensions=('time','y','x'), datatype='float32')

times.units = "hour"
xs.units = "m"
ys.units = "m"
rainfall.units = "mm"
times.axis = 'T'
xs.axis = 'X'
ys.axis = 'Y'

times[:] = np.linspace(0, (time_levels-1)*dt, time_levels)
xs[:] = blx + np.linspace(dx/2, xsz*dx - dx/2, xsz)
ys[:] = bly + np.linspace(dy/2, ysz*dy - dy/2, ysz)[::-1]

rainfall[0,:,:] = 0
rainfall[1,:,:] = 0
rainfall[1,1,1] = 10
rainfall[2,:,:] = 5

outfile.close()
