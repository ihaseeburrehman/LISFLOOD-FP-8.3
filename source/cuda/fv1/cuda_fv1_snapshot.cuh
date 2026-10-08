#pragma once
#include "cuda_snapshot.cuh"
#include "cuda_fv1_flow.cuh"

namespace lis
{
namespace cuda
{
namespace fv1
{

class Snapshot : public lis::cuda::Snapshot<Flow>
{
public:
	Snapshot
	(
		const char* resrootname,
		NUMERIC_TYPE interval,
		NUMERIC_TYPE next_save,
		int counter,
		NUMERIC_TYPE* DEM,
		Flow& U,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose,
		int precision = DEFAULT_PRECISION
	);

	void wait();

protected:
	void write(const int& call_gzip);
	void write_maxes();

private:
	const NUMERIC_TYPE* DEM;
	std::future<void> writer_elev;
	std::future<void> writer_H;
	std::future<void> writer_HU;
	std::future<void> writer_HV;
	std::future<void> writer_U;
	std::future<void> writer_V;
};

}
}
}
