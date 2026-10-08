#include "init_neighbours.cuh"
#include <helper_cuda.h> // MKS
#include <cub/cub.cuh> // MKS

__host__ lis::cuda::acc_nugrid::NonUniformNeighbours lis::cuda::acc_nugrid::init_neighbours
(
	AssembledSolution d_assem_sol,
	bool non_uniform_n
)
{
	NUMERIC_TYPE* h_num_neighbours = new NUMERIC_TYPE;
	NUMERIC_TYPE* d_num_neighbours = (NUMERIC_TYPE*)malloc_device(sizeof(NUMERIC_TYPE));

	void* d_temp_storage = NULL;
	size_t temp_storage = 0;

	CHECK_CUDA_ERROR(cub::DeviceReduce::Sum
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.nghbr_counts,
		d_num_neighbours,
		d_assem_sol.length
	));

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR(cub::DeviceReduce::Sum
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.nghbr_counts,
		d_num_neighbours,
		d_assem_sol.length
	));

	copy_cuda
	(
		h_num_neighbours, // dst
		d_num_neighbours, // src
		sizeof(NUMERIC_TYPE)
	);

	int num_neighbours = *h_num_neighbours;

	free_device(d_num_neighbours);
	free_device(d_temp_storage);
	delete      h_num_neighbours;

	// FINDING CUMULATIVE NGHBR COUNTS //

	d_temp_storage = NULL;
	temp_storage = 0;

	CHECK_CUDA_ERROR(cub::DeviceScan::ExclusiveSum
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.nghbr_counts,
		d_assem_sol.cumu_nghbr_counts,
		d_assem_sol.length
	));

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR(cub::DeviceScan::ExclusiveSum
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.nghbr_counts,
		d_assem_sol.cumu_nghbr_counts,
		d_assem_sol.length
	));

	free_device(d_temp_storage);

	// ------------------------------- //

	return NonUniformNeighbours(num_neighbours, non_uniform_n);
}