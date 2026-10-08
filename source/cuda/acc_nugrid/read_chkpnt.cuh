#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "AssembledSolution.h"
#include "NonUniformNeighbours.h"
#include "cuda_utils.cuh"
#include "compact.cuh"
#include "get_lvl_idx.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ void read_chkpnt
(
	AssembledSolution& d_assem_sol,
	NonUniformNeighbours& d_non_uniform_nghbrs,
	Fnames& filenames,
	States& states,
	Pars& pars,
	Solver& solver,
	NUMERIC_TYPE* d_maxH,
	NUMERIC_TYPE* d_totalHtm,
	NUMERIC_TYPE* d_maxHtm,
	NUMERIC_TYPE* d_initHtm,
	bool non_uniform_n,
	const int& num_finest_elems,
	int verbose
);

}
}
}