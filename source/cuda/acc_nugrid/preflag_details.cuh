#pragma once

#include "cuda_utils.cuh"
#include "get_lvl_idx.cuh"
#include "compact.cuh"
#include "Boundaries.h"
#include "PointSources.h"
#include "GaugePoints.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

bool* preflag_details
(
	const Boundaries&   boundaries,
	const PointSources& point_sources,
	const GaugePoints&  gauge_points,
	const int&          num_details,
	const int&          max_ref_lvl
);

}
}
}