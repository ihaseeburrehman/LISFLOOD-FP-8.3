#pragma once
#include "cuda_snapshot.cuh"
#include "cuda_dg2_dem.cuh"
#include "cuda_dg2_flow.cuh"

namespace lis
{
namespace cuda
{
namespace dg2
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
		Topography& DEM,
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
	const Topography& DEM;
	Flow& U;
	std::future<void> writer_elev;
	std::future<void> writer_H;
	std::future<void> writer_H1x;
	std::future<void> writer_H1y;
	std::future<void> writer_HU;
	std::future<void> writer_HU1x;
	std::future<void> writer_HU1y;
	std::future<void> writer_HV;
	std::future<void> writer_HV1x;
	std::future<void> writer_HV1y;
	std::future<void> writer_U;
	std::future<void> writer_V;
};

}
}
}
