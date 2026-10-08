
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "AssembledSolution.h"
#include "compact.cuh"
#include "cuda_utils.cuh"
#include "GaugePoints.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

void write_mass_data
(
	const char* respath,
//	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const Solver& solver,
	const States& states,
	const NUMERIC_TYPE& FArea,
	const NUMERIC_TYPE& vol2,
	const NUMERIC_TYPE& Qin,
	const NUMERIC_TYPE& Hds,
	const NUMERIC_TYPE& Qout,
	const NUMERIC_TYPE& Qerror,
	const NUMERIC_TYPE& Verror,
	const NUMERIC_TYPE& RainTotalLoss,
	const NUMERIC_TYPE& InfilTotalLoss,
	const NUMERIC_TYPE& EvapTotalLoss
);

}
}
}