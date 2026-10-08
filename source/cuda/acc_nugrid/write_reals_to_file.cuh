#pragma once

#include "cuda_runtime.h"

#include <stdio.h>
#include <string.h>

#include "cuda_utils.cuh"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void write_reals_to_file
(
	const char* filename,
	const char* respath,
	NUMERIC_TYPE*       d_results,
	const int&  array_length
);

}
}
}