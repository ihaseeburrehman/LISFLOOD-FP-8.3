
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


void write_gauge_point_data
(
	const char* respath,
	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const char* stagefilename,
	const Pars& pars,
	const Solver& solver,
	const States& states
);

void write_velocity_point_data
(
	const char* respath,
	AssembledSolution& d_assem_sol,
	const GaugePoints& gauge_points,
	const int& mesh_dim,
	const char* stagefilename,
	const Pars& pars,
	const Solver& solver,
	const States& states
);

}
}
}