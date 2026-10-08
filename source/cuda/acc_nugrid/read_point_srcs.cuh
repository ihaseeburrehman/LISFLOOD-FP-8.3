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

void read_point_srcs
(
	const char* bcifilename,
	const char* bdyfilename,
	const Pars& pars,
	const NUMERIC_TYPE& time_now,
	PointSources& point_sources,
	const BoundCs& h_src
);

}
}
}