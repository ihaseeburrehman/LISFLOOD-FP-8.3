#pragma once 

#include "cub/device/device_reduce.cuh"
#include "cuda_utils.cuh"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ NUMERIC_TYPE get_dt_CFL
(
	NUMERIC_TYPE*&     d_dt_CFL, 
	const int& sol_len
);

}
}
}