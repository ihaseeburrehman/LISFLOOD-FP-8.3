#pragma once

#include "cuda_runtime.h"
#include "../../lisflood.h"
#include "ScaleChildren.h"
#include "maskFilterHaar.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

//__device__ __forceinline__ NUMERIC_TYPE encode_scale(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * (H0 * (H0 * u.child_0 + H1 * u.child_2) + H1 * (H0 * u.child_1 + H1 * u.child_3));
//}

/*__device__ __forceinline__ NUMERIC_TYPE encode_scale_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
{
	return C(0.5) * (H0_11 * (H0_11 * u0.child_0  + H0_12 * u1y.child_0 + H1_11 * u0.child_2  + H1_12 * u1y.child_2) +
		             H0_12 * (H0_11 * u1x.child_0 + H0_12 * uxy.child_0 + H1_11 * u1x.child_2 + H1_12 * uxy.child_2) +
		             H1_11 * (H0_11 * u0.child_1  + H0_12 * u1y.child_1 + H1_11 * u0.child_3  + H1_12 * u1y.child_3) +
		             H1_12 * (H0_11 * u1x.child_1 + H0_12 * uxy.child_1 + H1_11 * u1x.child_3 + H1_12 * uxy.child_3));
}

__device__ __forceinline__ NUMERIC_TYPE encode_scale_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
{
	return C(0.5) * (H0_11 * (H0_21 * u0.child_0  + H0_22 * u1y.child_0 + H1_21 * u0.child_2  + H1_22 * u1y.child_2) +
		             H0_12 * (H0_21 * u1x.child_0 + H0_22 * uxy.child_0 + H1_21 * u1x.child_2 + H1_22 * uxy.child_2) +
		             H1_11 * (H0_21 * u0.child_1  + H0_22 * u1y.child_1 + H1_21 * u0.child_3  + H1_22 * u1y.child_3) +
		             H1_12 * (H0_21 * u1x.child_1 + H0_22 * uxy.child_1 + H1_21 * u1x.child_3 + H1_22 * uxy.child_3));
}

__device__ __forceinline__ NUMERIC_TYPE encode_scale_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
{
	return C(0.5) * (H0_21 * (H0_11 * u0.child_0  + H0_12 * u1y.child_0 + H1_11 * u0.child_2  + H1_12 * u1y.child_2) +
		             H0_22 * (H0_11 * u1x.child_0 + H0_12 * uxy.child_0 + H1_11 * u1x.child_2 + H1_12 * uxy.child_2) +
		             H1_21 * (H0_11 * u0.child_1  + H0_12 * u1y.child_1 + H1_11 * u0.child_3  + H1_12 * u1y.child_3) +
		             H1_22 * (H0_11 * u1x.child_1 + H0_12 * uxy.child_1 + H1_11 * u1x.child_3 + H1_12 * uxy.child_3));
}

__device__ __forceinline__ NUMERIC_TYPE encode_scale_22(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
{
	return C(0.5) * (H0_21 * (H0_21 * u0.child_0  + H0_22 * u1y.child_0 + H1_21 * u0.child_2  + H1_22 * u1y.child_2) +
		             H0_22 * (H0_21 * u1x.child_0 + H0_22 * uxy.child_0 + H1_21 * u1x.child_2 + H1_22 * uxy.child_2) +
		             H1_21 * (H0_21 * u0.child_1  + H0_22 * u1y.child_1 + H1_21 * u0.child_3  + H1_22 * u1y.child_3) +
		             H1_22 * (H0_21 * u1x.child_1 + H0_22 * uxy.child_1 + H1_21 * u1x.child_3 + H1_22 * uxy.child_3));
}*/		

__device__ __forceinline__ NUMERIC_TYPE encode_scale_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return HH0_11 * u0.child_0 + HH0_12 * u1x.child_0 + HH0_13 * u1y.child_0 +
		   HH1_11 * u0.child_2 + HH1_12 * u1x.child_2 + HH1_13 * u1y.child_2 +
		   HH2_11 * u0.child_1 + HH2_12 * u1x.child_1 + HH2_13 * u1y.child_1 +
		   HH3_11 * u0.child_3 + HH3_12 * u1x.child_3 + HH3_13 * u1y.child_3;
}

__device__ __forceinline__ NUMERIC_TYPE encode_scale_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return HH0_21 * u0.child_0 + HH0_22 * u1x.child_0 + HH0_23 * u1y.child_0 +
		   HH1_21 * u0.child_2 + HH1_22 * u1x.child_2 + HH1_23 * u1y.child_2 +
		   HH2_21 * u0.child_1 + HH2_22 * u1x.child_1 + HH2_23 * u1y.child_1 +
		   HH3_21 * u0.child_3 + HH3_22 * u1x.child_3 + HH3_23 * u1y.child_3;
}

__device__ __forceinline__ NUMERIC_TYPE encode_scale_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return HH0_31 * u0.child_0 + HH0_32 * u1x.child_0 + HH0_33 * u1y.child_0 +
		   HH1_31 * u0.child_2 + HH1_32 * u1x.child_2 + HH1_33 * u1y.child_2 +
		   HH2_31 * u0.child_1 + HH2_32 * u1x.child_1 + HH2_33 * u1y.child_1 +
		   HH3_31 * u0.child_3 + HH3_32 * u1x.child_3 + HH3_33 * u1y.child_3;
}


}
}
}