#pragma once
#include "../../lisflood.h"
#include "../../geometry.h"
#include "../cuda_simulate.cuh"

namespace lis
{
namespace cuda
{
namespace adaptive
{

class Simulation : public cuda::Simulation
{
public:
	void run
	(
		int    argc,
		char** argv
	);
};

}
}
}

