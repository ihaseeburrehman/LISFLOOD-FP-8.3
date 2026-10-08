#pragma once

#pragma once
#include "../../lisflood.h"
#include "../../geometry.h"
#include "../cuda_simulate.cuh"

namespace lis
{
namespace cuda
{
namespace dg2
{

class Simulation : public cuda::Simulation
{
public:
	void run
	(
		Fnames& filenames,
		States& states,
		Pars& pars,
		::Solver& solver,
		int verbose
	);

private:
	void initialise_H
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* H1x,
		NUMERIC_TYPE* H1y,
		const char* filename,
		States& states,
		NUMERIC_TYPE* DEM,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void initialise_H_slope
	(
		NUMERIC_TYPE* array,
		const char* filename,
		const char* suffix,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);
};

}
}
}

