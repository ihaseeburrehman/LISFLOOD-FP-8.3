#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "generate_morton_code.cuh"
#include "../../lisflood.h"
#include "GaugePoints.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

void read_gauge_points
(
	const char* stagefilename,
	const Pars& pars,
	const GaugePoints& gauge_points
);

}
}
}