#include "ghostraster.h"
#include "cuda_geometry.cuh"
#include "cuda_util.cuh"
#include <helper_cuda.h>
#include <cub/cub.cuh>

template<typename F>
lis::cuda::DynamicTimestep<F>::DynamicTimestep
(
	NUMERIC_TYPE& dt,
	Geometry& geometry,
	NUMERIC_TYPE max_dt,
	int adaptive_ts,
	Solver<F>& solver
)
:
dt(dt),
max_dt(max_dt),
adaptive(adaptive_ts == ON),
solver(solver),
elements(lis::GhostRaster::elements(geometry))
{
	if (adaptive)
	{
		d_temp = nullptr;
		NUMERIC_TYPE* dummy_in = nullptr;
		NUMERIC_TYPE* dummy_out = nullptr;

		checkCudaErrors(cub::DeviceReduce::Min(d_temp, bytes,
					dummy_in, dummy_out, elements));

		d_temp = malloc_device(bytes);

		dt_field = cuda::GhostRaster::allocate_device(geometry);
	}
	else
	{
		dt = max_dt;
	}
}

template<typename F>
lis::cuda::DynamicTimestepACC<F>::DynamicTimestepACC
(
	NUMERIC_TYPE& dt,
	Geometry& geometry,
	NUMERIC_TYPE max_dt,
	int adaptive_ts,
	Solver<F>& solver
)
	:
	dt(dt),
	max_dt(max_dt),
	adaptive(adaptive_ts == ON),
	solver(solver),
	elements(lis::GhostRaster::elements_H(geometry))
{

		d_temp = nullptr;
		NUMERIC_TYPE* dummy_in = nullptr;
		NUMERIC_TYPE* dummy_out = nullptr;

		checkCudaErrors(cub::DeviceReduce::Min(d_temp, bytes,
			dummy_in, dummy_out, elements));

		d_temp = malloc_device(bytes);

		dt_field = cuda::GhostRaster::allocate_device_H(geometry);


}

template<typename F>
NUMERIC_TYPE lis::cuda::DynamicTimestep<F>::update_dt()
{
	if (adaptive)
	{
		solver.update_dt_per_element(dt_field);
		checkCudaErrors(cub::DeviceReduce::Min(d_temp, bytes,
					dt_field, &dt, elements));
	}
	cuda::sync();

	return dt;
}

template<typename F>
NUMERIC_TYPE lis::cuda::DynamicTimestepACC<F>::update_dt_ACC()
{
		solver.update_dt_per_element(dt_field);
		checkCudaErrors(cub::DeviceReduce::Min(d_temp, bytes,
					dt_field, &dt, elements));

	cuda::sync(); 

	return dt;
}

template<typename F>
lis::cuda::DynamicTimestep<F>::~DynamicTimestep()
{
	if (adaptive)
	{
		cuda::free_device(dt_field);
		free_device(d_temp);
	}
}

template<typename F>
lis::cuda::DynamicTimestepACC<F>::~DynamicTimestepACC()
{

		cuda::free_device(dt_field);
		free_device(d_temp);

}