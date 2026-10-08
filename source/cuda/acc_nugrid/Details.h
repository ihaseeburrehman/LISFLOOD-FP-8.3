#pragma once

#include "SubDetails.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct Details
{
//	SubDetails eta0;
//	SubDetails qx0;
//	SubDetails qy0;
	SubDetails z0;

//	SubDetails eta1x;
//	SubDetails qx1x;
//	SubDetails qy1x;
	SubDetails z1x;
	
//	SubDetails eta1y;
//	SubDetails qx1y;
//	SubDetails qy1y;
	SubDetails z1y;
//	SubDetails zxy;

	SubDetails n0;
	SubDetails n1x;
	SubDetails n1y;

	SubDetails h0;
	SubDetails h1x;
	SubDetails h1y;

	Details
	(
		const int& num_details, 
		bool non_uniform_n,
		const int& startfile
	) 
	: 
		z0(num_details), 
		z1x(num_details), 
		z1y(num_details), 
		n0(non_uniform_n ? num_details : 0),
		n1x(non_uniform_n ? num_details : 0),
		n1y(non_uniform_n ? num_details : 0),
		h0(startfile ? num_details : 0),
		h1x(startfile ? num_details : 0),
		h1y(startfile ? num_details : 0)
	{}

} Details;

}
}
}