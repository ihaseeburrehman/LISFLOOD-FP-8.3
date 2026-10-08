#pragma once

#include "cuda_runtime.h"
#include <math.h>
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct SubDetail
{
	NUMERIC_TYPE alpha_11;
	NUMERIC_TYPE alpha_21;
	NUMERIC_TYPE alpha_12;

	NUMERIC_TYPE beta_11;
	NUMERIC_TYPE beta_21;
	NUMERIC_TYPE beta_12;

	NUMERIC_TYPE gamma_11;
	NUMERIC_TYPE gamma_21;
	NUMERIC_TYPE gamma_12;

	__device__  __forceinline__ NUMERIC_TYPE get_max()
	{
		//NUMERIC_TYPE max_detail = FMAX( FABS(alpha), FABS(beta) );
		//     max_detail = FMAX( FABS(gamma), max_detail);

		NUMERIC_TYPE max_detail = FMAX(FABS(alpha_11), FABS(alpha_12));
		max_detail = FMAX(FABS(alpha_21), max_detail);
		max_detail = FMAX(FABS(beta_11), max_detail);
		max_detail = FMAX(FABS(beta_12), max_detail);
		max_detail = FMAX(FABS(beta_21), max_detail);
		max_detail = FMAX(FABS(gamma_11), max_detail);
		max_detail = FMAX(FABS(gamma_12), max_detail);
		max_detail = FMAX(FABS(gamma_21), max_detail);

		return max_detail;
	}
} SubDetail;

}
}
}