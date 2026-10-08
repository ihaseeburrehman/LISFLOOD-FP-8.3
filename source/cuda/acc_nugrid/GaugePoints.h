#pragma once

#include "MortonCode.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct GaugePoints
{
	MortonCode* codes;
	int         num_points;
	bool        is_copy = false;

	GaugePoints(const int& num_points)
	{
		this->codes = (num_points > 0) ? new MortonCode[num_points] : nullptr;
		this->num_points = num_points;
	}

	GaugePoints(const GaugePoints& original) { *this = original; is_copy = true; }

	~GaugePoints() { if (!is_copy && codes != nullptr) delete[] codes; }

} GaugePoints;

}
}
}