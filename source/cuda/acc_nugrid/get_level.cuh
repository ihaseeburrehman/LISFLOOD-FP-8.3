#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include "index_1D.h"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__device__ int get_level(index_1D idx);

}
}
}