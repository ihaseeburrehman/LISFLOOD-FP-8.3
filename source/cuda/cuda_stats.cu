#include "cuda_stats.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include <cub/cub.cuh>
#include <limits>

namespace lis
{
namespace cuda
{

__global__ void accumulate_mass
(
	MassStats* cumulative,
	MassStats* instantaneous
)
{
	cumulative->in += instantaneous->in * cuda::dt;
	cumulative->out += instantaneous->out * cuda::dt;
}

}
}

lis::cuda::StatsCollector::StatsCollector
(
	Geometry& geometry,
	NUMERIC_TYPE DepthThresh,
	int acceleration
)
:
geometry(geometry),
elements((acceleration == 0) ? lis::GhostRaster::elements(geometry) : lis::GhostRaster::elements_H(geometry))
{
	area_temp = nullptr;
	volume_temp = nullptr;
	NUMERIC_TYPE* dummy_in = nullptr;
	NUMERIC_TYPE* dummy_reduce_out = nullptr;
	int* dummy_histogram_out = nullptr;

	histogram_levels = static_cast<NUMERIC_TYPE*>(
			malloc_unified(3*sizeof(NUMERIC_TYPE)));

	histogram_levels[0] = C(0.0);
	histogram_levels[1] = DepthThresh;
	histogram_levels[2] = std::numeric_limits<NUMERIC_TYPE>::max();

	histogram_counts = static_cast<int*>(malloc_unified(2*sizeof(int)));

	checkCudaErrors(cub::DeviceHistogram::HistogramRange(area_temp, area_bytes,
				dummy_in, dummy_histogram_out, 3, histogram_levels, elements));
	area_temp = malloc_device(area_bytes);

	checkCudaErrors(cub::DeviceReduce::Sum(volume_temp, volume_bytes,
				dummy_in, dummy_reduce_out, elements));
	volume_temp = malloc_device(volume_bytes);

	d_volume = static_cast<NUMERIC_TYPE*>(malloc_device(sizeof(NUMERIC_TYPE)));

	d_instantaneous_mass = static_cast<MassStats*>(
			malloc_device(sizeof(MassStats)));
	d_cumulative_mass = static_cast<MassStats*>(
			malloc_device(sizeof(MassStats)));
	zero_cumulative_mass();
}

void lis::cuda::StatsCollector::zero_instantaneous_mass()
{
	MassStats temp = {};
	cuda::copy(d_instantaneous_mass, &temp, sizeof(MassStats));
}

lis::MassStats* lis::cuda::StatsCollector::instantaneous_mass()
{
	return d_instantaneous_mass;
}

void lis::cuda::StatsCollector::accumulate_mass()
{
	lis::cuda::accumulate_mass<<<1, 1>>>(d_cumulative_mass,
			d_instantaneous_mass);
	cuda::sync();
	cumulative_time += cuda::dt;
}

lis::StatsEntry lis::cuda::StatsCollector::create_entry
(
	::Solver& solver,
	NUMERIC_TYPE* H
)
{
	NUMERIC_TYPE current_volume = volume(H);
	NUMERIC_TYPE volume_change = current_volume - previous_volume;
	previous_volume = current_volume;

	MassStats h_instantaneous_mass;
	cuda::copy(&h_instantaneous_mass, d_instantaneous_mass, sizeof(MassStats));

	MassStats h_cumulative_mass;
	cuda::copy(&h_cumulative_mass, d_cumulative_mass, sizeof(MassStats));

	NUMERIC_TYPE time = cumulative_time;

	NUMERIC_TYPE net_source_flux =
		h_cumulative_mass.in - h_cumulative_mass.out;

	zero_cumulative_mass();

	return { 
		solver.t,
		cuda::dt,
		solver.itCount,
		solver.MinTstep,
		area(H),
		current_volume,
		h_instantaneous_mass,
		net_source_flux - volume_change,
		(net_source_flux - volume_change)/time
	};
}

lis::cuda::StatsCollector::~StatsCollector()
{
	free_unified(histogram_levels);
	free_unified(histogram_counts);
	free_device(area_temp);
	free_device(volume_temp);
	free_device(d_volume);
	free_device(d_instantaneous_mass);
	free_device(d_cumulative_mass);
}

NUMERIC_TYPE lis::cuda::StatsCollector::area
(
	NUMERIC_TYPE* H
)
{
	checkCudaErrors(cub::DeviceHistogram::HistogramRange(area_temp, area_bytes,
				H, histogram_counts, 3, histogram_levels, elements));
	cuda::sync();
	return histogram_counts[1] * geometry.dx * geometry.dy;
}

NUMERIC_TYPE lis::cuda::StatsCollector::volume
(
	NUMERIC_TYPE* H
)
{

	checkCudaErrors(cub::DeviceReduce::Sum(volume_temp, volume_bytes,
				H, d_volume, elements));

	NUMERIC_TYPE h_volume;

	cuda::copy(&h_volume, d_volume, sizeof(NUMERIC_TYPE));

	return h_volume * geometry.dx * geometry.dy;
}

void lis::cuda::StatsCollector::zero_cumulative_mass()
{
	MassStats temp_stats = {};
	cuda::copy(d_cumulative_mass, &temp_stats, sizeof(MassStats));

	cumulative_time = C(0.0);
}
