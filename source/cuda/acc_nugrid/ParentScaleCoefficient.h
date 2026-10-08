#pragma once

#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct ParentScaleCoefficient
{
//	NUMERIC_TYPE eta;
//	NUMERIC_TYPE qx;
//	NUMERIC_TYPE qy;
	NUMERIC_TYPE z0;
	NUMERIC_TYPE z1x;
	NUMERIC_TYPE z1y;
	NUMERIC_TYPE zxy;

} ParentScaleCoefficient;

}
}
}