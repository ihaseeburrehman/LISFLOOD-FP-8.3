#pragma once
#include "../lisflood.h"
#include "../geometry.h"

namespace lis
{
namespace cuda
{

class MaxField
{
public:
	MaxField
	(
		const char* resrootname,
		Geometry& geometry,
		int pitch,
		int offset,
		dim3 grid_size,
		int acceleration,
		int precision = DEFAULT_PRECISION
		
	);

	void update(NUMERIC_TYPE* H);
	void updateACC(NUMERIC_TYPE* H);
	void write();

	~MaxField();

private:
	const char* resrootname;
	Geometry& geometry;
	int pitch;
	int offset;
	dim3 grid_size;
	int precision;
	NUMERIC_TYPE* maxH;
	NUMERIC_TYPE* totalHtm;
	NUMERIC_TYPE* maxHtm;
	NUMERIC_TYPE* initHtm;
};

}
}
