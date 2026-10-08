#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "AssembledSolution.h"
#include "cuda_utils.cuh"
#include "compact.cuh"
#include "get_lvl_idx.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct Points
{
	NUMERIC_TYPE ll_x;
	NUMERIC_TYPE ll_y;
	NUMERIC_TYPE ul_x;
	NUMERIC_TYPE ul_y;
	NUMERIC_TYPE lr_x;
	NUMERIC_TYPE lr_y;
	NUMERIC_TYPE ur_x;
	NUMERIC_TYPE ur_y;
	
} Points;

__host__ void write_soln_vtk
(
	const char*                 respath,
	const AssembledSolution&    d_assem_sol,
	const Pars& pars,
	const int&                  lev,
	const NUMERIC_TYPE&                 depth_thresh,
	const int& call_gzip,
	const int& precision
);

__host__ void write_soln_vtk_with_n
(
	const char* respath,
	const AssembledSolution& d_assem_sol,
	const Pars& pars,
	const int& lev,
	const NUMERIC_TYPE& depth_thresh,
	const int& call_gzip,
	const int& precision
);

}
}
}