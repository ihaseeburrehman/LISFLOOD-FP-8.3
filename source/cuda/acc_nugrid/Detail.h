#pragma once

#include "SubDetail.h"
#include "Maxes.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct Detail
{
//	SubDetail eta;
//	SubDetail qx;
//	SubDetail qy;
	SubDetail z;

	__device__ __forceinline__
	NUMERIC_TYPE get_normalised_detail(Maxes maxes)
	{
		NUMERIC_TYPE normalised_detail = C(0.0);
		NUMERIC_TYPE z_normalised   = z.get_max()   / maxes.z;
		normalised_detail = z_normalised;
		return normalised_detail;
	}

} Detail;

}
}
}