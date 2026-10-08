#pragma once

#include "cuda_utils.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PointSources.h"
#include "generate_morton_code.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

PointSources read_num_point_srcs
(
	const char* bcifilename
);

}
}
}