#pragma once

#pragma once
#include "../../lisflood.h"
#include "../../geometry.h"
#include "../cuda_simulate.cuh"

namespace lis
{
namespace cuda
{
namespace fv2
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
};

}
}
}

