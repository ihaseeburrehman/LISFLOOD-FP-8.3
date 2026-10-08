#pragma once
#include "../../lisflood.h"
#include "../../geometry.h"
#include "../cuda_simulate.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
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

