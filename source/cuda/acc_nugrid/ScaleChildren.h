#pragma once

#include "cuda_runtime.h"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct ScaleChildren
{
	NUMERIC_TYPE child_0;
	NUMERIC_TYPE child_1;
	NUMERIC_TYPE child_2;
	NUMERIC_TYPE child_3;

} ScaleChildren;

}
}
}