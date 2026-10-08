#pragma once

#include "cuda_runtime.h"
#include "Detail.h"
#include "ChildScaleCoefficients.h"
#include "maskFilterHaar.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

// encodes the details alpha, beta and gamma for eta, qx, qy and z
__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);

__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);

__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);
__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy);

//__device__ __forceinline__ Detail encode_details(ChildScaleCoefficients& child_scale_coeffs)
//{
//	SubDetail eta =
//	{
//		encode_detail_alpha(child_scale_coeffs.eta),
//		encode_detail_beta(child_scale_coeffs.eta),
//		encode_detail_gamma(child_scale_coeffs.eta)
//	};
//
//	SubDetail qx =
//	{
//		encode_detail_alpha(child_scale_coeffs.qx),
//		encode_detail_beta(child_scale_coeffs.qx),
//		encode_detail_gamma(child_scale_coeffs.qx)
//	};
//
//	SubDetail qy =
//	{
//		encode_detail_alpha(child_scale_coeffs.qy),
//		encode_detail_beta(child_scale_coeffs.qy),
//		encode_detail_gamma(child_scale_coeffs.qy)
//	};
//	
//	SubDetail z =
//	{
//		encode_detail_alpha(child_scale_coeffs.z),
//		encode_detail_beta(child_scale_coeffs.z),
//		encode_detail_gamma(child_scale_coeffs.z)
//	};
//
//	Detail details =
//	{
//		eta,
//		qx,
//		qy,
//		z
//	};
//
//	return details;
//}

//__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( H0_11 * (G0_11 * u0.child_0  + G0_12 * u1y.child_0 + G1_11 * u0.child_2  + G1_12 * u1y.child_2 ) +
//		              H0_12 * (G0_11 * u1x.child_0 + G0_12 * uxy.child_0 + G1_11 * u1x.child_2 + G1_12 * uxy.child_2 ) + 
//		              H1_11 * (G0_11 * u0.child_1  + G0_12 * u1y.child_1 + G1_11 * u0.child_3  + G1_12 * u1y.child_3 ) + 
//		              H1_12 * (G0_11 * u1x.child_1 + G0_12 * uxy.child_1 + G1_11 * u1x.child_3 + G1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( H0_11 * (G0_21 * u0.child_0  + G0_22 * u1y.child_0 + G1_21 * u0.child_2  + G1_22 * u1y.child_2 ) +
//		              H0_12 * (G0_21 * u1x.child_0 + G0_22 * uxy.child_0 + G1_21 * u1x.child_2 + G1_22 * uxy.child_2 ) + 
//		              H1_11 * (G0_21 * u0.child_1  + G0_22 * u1y.child_1 + G1_21 * u0.child_3  + G1_22 * u1y.child_3 ) + 
//		              H1_12 * (G0_21 * u1x.child_1 + G0_22 * uxy.child_1 + G1_21 * u1x.child_3 + G1_22 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( H0_21 * (G0_11 * u0.child_0  + G0_12 * u1y.child_0 + G1_11 * u0.child_2  + G1_12 * u1y.child_2 ) +
//		              H0_22 * (G0_11 * u1x.child_0 + G0_12 * uxy.child_0 + G1_11 * u1x.child_2 + G1_12 * uxy.child_2 ) + 
//		              H1_21 * (G0_11 * u0.child_1  + G0_12 * u1y.child_1 + G1_11 * u0.child_3  + G1_12 * u1y.child_3 ) + 
//		              H1_22 * (G0_11 * u1x.child_1 + G0_12 * uxy.child_1 + G1_11 * u1x.child_3 + G1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_22(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( H0_21 * (G0_21 * u0.child_0  + G0_22 * u1y.child_0 + G1_21 * u0.child_2  + G1_22 * u1y.child_2 ) +
//		              H0_22 * (G0_21 * u1x.child_0 + G0_22 * uxy.child_0 + G1_21 * u1x.child_2 + G1_22 * uxy.child_2 ) + 
//		              H1_21 * (G0_21 * u0.child_1  + G0_22 * u1y.child_1 + G1_21 * u0.child_3  + G1_22 * u1y.child_3 ) + 
//		              H1_22 * (G0_21 * u1x.child_1 + G0_22 * uxy.child_1 + G1_21 * u1x.child_3 + G1_22 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_11 * (H0_11 * u0.child_0  + H0_12 * u1y.child_0 + H1_11 * u0.child_2  + H1_12 * u1y.child_2 ) +
//		              G0_12 * (H0_11 * u1x.child_0 + H0_12 * uxy.child_0 + H1_11 * u1x.child_2 + H1_12 * uxy.child_2 ) + 
//		              G1_11 * (H0_11 * u0.child_1  + H0_12 * u1y.child_1 + H1_11 * u0.child_3  + H1_12 * u1y.child_3 ) + 
//		              G1_12 * (H0_11 * u1x.child_1 + H0_12 * uxy.child_1 + H1_11 * u1x.child_3 + H1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_11 * (H0_21 * u0.child_0  + H0_22 * u1y.child_0 + H1_21 * u0.child_2  + H1_22 * u1y.child_2 ) +
//		              G0_12 * (H0_21 * u1x.child_0 + H0_22 * uxy.child_0 + H1_21 * u1x.child_2 + H1_22 * uxy.child_2 ) + 
//		              G1_11 * (H0_21 * u0.child_1  + H0_22 * u1y.child_1 + H1_21 * u0.child_3  + H1_22 * u1y.child_3 ) + 
//		              G1_12 * (H0_21 * u1x.child_1 + H0_22 * uxy.child_1 + H1_21 * u1x.child_3 + H1_22 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_21 * (H0_11 * u0.child_0  + H0_12 * u1y.child_0 + H1_11 * u0.child_2  + H1_12 * u1y.child_2 ) +
//		              G0_22 * (H0_11 * u1x.child_0 + H0_12 * uxy.child_0 + H1_11 * u1x.child_2 + H1_12 * uxy.child_2 ) + 
//		              G1_21 * (H0_11 * u0.child_1  + H0_12 * u1y.child_1 + H1_11 * u0.child_3  + H1_12 * u1y.child_3 ) + 
//		              G1_22 * (H0_11 * u1x.child_1 + H0_12 * uxy.child_1 + H1_11 * u1x.child_3 + H1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_22(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_21 * (H0_21 * u0.child_0  + H0_22 * u1y.child_0 + H1_21 * u0.child_2  + H1_22 * u1y.child_2 ) +
//		              G0_22 * (H0_21 * u1x.child_0 + H0_22 * uxy.child_0 + H1_21 * u1x.child_2 + H1_22 * uxy.child_2 ) + 
//		              G1_21 * (H0_21 * u0.child_1  + H0_22 * u1y.child_1 + H1_21 * u0.child_3  + H1_22 * u1y.child_3 ) + 
//		              G1_22 * (H0_21 * u1x.child_1 + H0_22 * uxy.child_1 + H1_21 * u1x.child_3 + H1_22 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_11 * (G0_11 * u0.child_0  + G0_12 * u1y.child_0 + G1_11 * u0.child_2  + G1_12 * u1y.child_2 ) +
//		              G0_12 * (G0_11 * u1x.child_0 + G0_12 * uxy.child_0 + G1_11 * u1x.child_2 + G1_12 * uxy.child_2 ) + 
//		              G1_11 * (G0_11 * u0.child_1  + G0_12 * u1y.child_1 + G1_11 * u0.child_3  + G1_12 * u1y.child_3 ) + 
//		              G1_12 * (G0_11 * u1x.child_1 + G0_12 * uxy.child_1 + G1_11 * u1x.child_3 + G1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_11 * (G0_21 * u0.child_0  + G0_22 * u1y.child_0 + G1_21 * u0.child_2  + G1_22 * u1y.child_2 ) +
//		              G0_12 * (G0_21 * u1x.child_0 + G0_22 * uxy.child_0 + G1_21 * u1x.child_2 + G1_22 * uxy.child_2 ) + 
//		              G1_11 * (G0_21 * u0.child_1  + G0_22 * u1y.child_1 + G1_21 * u0.child_3  + G1_22 * u1y.child_3 ) + 
//		              G1_12 * (G0_21 * u1x.child_1 + G0_22 * uxy.child_1 + G1_21 * u1x.child_3 + G1_22 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_21 * (G0_11 * u0.child_0  + G0_12 * u1y.child_0 + G1_11 * u0.child_2  + G1_12 * u1y.child_2 ) +
//		              G0_22 * (G0_11 * u1x.child_0 + G0_12 * uxy.child_0 + G1_11 * u1x.child_2 + G1_12 * uxy.child_2 ) + 
//		              G1_21 * (G0_11 * u0.child_1  + G0_12 * u1y.child_1 + G1_11 * u0.child_3  + G1_12 * u1y.child_3 ) + 
//		              G1_22 * (G0_11 * u1x.child_1 + G0_12 * uxy.child_1 + G1_11 * u1x.child_3 + G1_12 * uxy.child_3 ) );
//}
//
//__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_22(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y, ScaleChildren& uxy)
//{
//	return C(0.5) * ( G0_21 * (G0_21 * u0.child_0  + G0_22 * u1y.child_0 + G1_21 * u0.child_2  + G1_22 * u1y.child_2 ) +
//		              G0_22 * (G0_21 * u1x.child_0 + G0_22 * uxy.child_0 + G1_21 * u1x.child_2 + G1_22 * uxy.child_2 ) + 
//		              G1_21 * (G0_21 * u0.child_1  + G0_22 * u1y.child_1 + G1_21 * u0.child_3  + G1_22 * u1y.child_3 ) + 
//		              G1_22 * (G0_21 * u1x.child_1 + G0_22 * uxy.child_1 + G1_21 * u1x.child_3 + G1_22 * uxy.child_3 ) );
//}



__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GA0_11 * u0.child_0 + GA0_12 * u1x.child_0 + GA0_13 * u1y.child_0 +
		   GA1_11 * u0.child_2 + GA1_12 * u1x.child_2 + GA1_13 * u1y.child_2 +
		   GA2_11 * u0.child_1 + GA2_12 * u1x.child_1 + GA2_13 * u1y.child_1 +
		   GA3_11 * u0.child_3 + GA3_12 * u1x.child_3 + GA3_13 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GA0_21 * u0.child_0 + GA0_22 * u1x.child_0 + GA0_23 * u1y.child_0 +
		   GA1_21 * u0.child_2 + GA1_22 * u1x.child_2 + GA1_23 * u1y.child_2 +
		   GA2_21 * u0.child_1 + GA2_22 * u1x.child_1 + GA2_23 * u1y.child_1 +
		   GA3_21 * u0.child_3 + GA3_22 * u1x.child_3 + GA3_23 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_alpha_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GA0_31 * u0.child_0 + GA0_32 * u1x.child_0 + GA0_33 * u1y.child_0 +
		   GA1_31 * u0.child_2 + GA1_32 * u1x.child_2 + GA1_33 * u1y.child_2 +
		   GA2_31 * u0.child_1 + GA2_32 * u1x.child_1 + GA2_33 * u1y.child_1 +
		   GA3_31 * u0.child_3 + GA3_32 * u1x.child_3 + GA3_33 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GB0_11 * u0.child_0 + GB0_12 * u1x.child_0 + GB0_13 * u1y.child_0 +
		   GB1_11 * u0.child_2 + GB1_12 * u1x.child_2 + GB1_13 * u1y.child_2 +
		   GB2_11 * u0.child_1 + GB2_12 * u1x.child_1 + GB2_13 * u1y.child_1 +
		   GB3_11 * u0.child_3 + GB3_12 * u1x.child_3 + GB3_13 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GB0_21 * u0.child_0 + GB0_22 * u1x.child_0 + GB0_23 * u1y.child_0 +
		   GB1_21 * u0.child_2 + GB1_22 * u1x.child_2 + GB1_23 * u1y.child_2 +
		   GB2_21 * u0.child_1 + GB2_22 * u1x.child_1 + GB2_23 * u1y.child_1 +
		   GB3_21 * u0.child_3 + GB3_22 * u1x.child_3 + GB3_23 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_beta_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GB0_31 * u0.child_0 + GB0_32 * u1x.child_0 + GB0_33 * u1y.child_0 +
		   GB1_31 * u0.child_2 + GB1_32 * u1x.child_2 + GB1_33 * u1y.child_2 +
		   GB2_31 * u0.child_1 + GB2_32 * u1x.child_1 + GB2_33 * u1y.child_1 +
		   GB3_31 * u0.child_3 + GB3_32 * u1x.child_3 + GB3_33 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_11(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GC0_11 * u0.child_0 + GC0_12 * u1x.child_0 + GC0_13 * u1y.child_0 +
		   GC1_11 * u0.child_2 + GC1_12 * u1x.child_2 + GC1_13 * u1y.child_2 +
		   GC2_11 * u0.child_1 + GC2_12 * u1x.child_1 + GC2_13 * u1y.child_1 +
		   GC3_11 * u0.child_3 + GC3_12 * u1x.child_3 + GC3_13 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_21(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GC0_21 * u0.child_0 + GC0_22 * u1x.child_0 + GC0_23 * u1y.child_0 +
		   GC1_21 * u0.child_2 + GC1_22 * u1x.child_2 + GC1_23 * u1y.child_2 +
		   GC2_21 * u0.child_1 + GC2_22 * u1x.child_1 + GC2_23 * u1y.child_1 +
		   GC3_21 * u0.child_3 + GC3_22 * u1x.child_3 + GC3_23 * u1y.child_3 ;
}

__device__ __forceinline__ NUMERIC_TYPE encode_detail_gamma_12(ScaleChildren& u0, ScaleChildren& u1x, ScaleChildren& u1y)
{
	return GC0_31 * u0.child_0 + GC0_32 * u1x.child_0 + GC0_33 * u1y.child_0 +
		   GC1_31 * u0.child_2 + GC1_32 * u1x.child_2 + GC1_33 * u1y.child_2 +
		   GC2_31 * u0.child_1 + GC2_32 * u1x.child_1 + GC2_33 * u1y.child_1 +
		   GC3_31 * u0.child_3 + GC3_32 * u1x.child_3 + GC3_33 * u1y.child_3 ;
}

}
}
}