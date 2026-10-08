#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

void read_raster_file
(
	const char* raster_filename,
	NUMERIC_TYPE* raster_array,
	const int&  extended_mesh_dim,
	const NUMERIC_TYPE& no_data
);

}
}
}