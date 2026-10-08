#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_bdy_conds.cuh"
#include "Boundaries.h"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

	void read_all_bdy_conds
(
	const char* bcifilename,
	Boundaries& boundaries,
	const Pars& pars
);

}
}
}