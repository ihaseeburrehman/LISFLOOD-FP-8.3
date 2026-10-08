#pragma once

#include "ScaleChildren.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct ChildScaleCoefficients
{
//	ScaleChildren eta;
//	ScaleChildren qx;
//	ScaleChildren qy;
	ScaleChildren z0;
	ScaleChildren z1x;
	ScaleChildren z1y;
	ScaleChildren zxy;

} ChildScaleCoefficients;

}
}
}