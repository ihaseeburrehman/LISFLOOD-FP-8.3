#include "sort_finest_scale_coefficients_z_order.cuh"

__host__ void lis::cuda::acc_nugrid::sort_finest_scale_coefficients_z_order
(
	MortonCode*        d_morton_codes,
	MortonCode*        d_sorted_morton_codes,
	AssembledSolution& d_assembled_solution,
	AssembledSolution& d_buffer_assembled_solution,
	MortonCode*        d_indices,
	MortonCode*        d_reverse_z_order,
	bool               non_uniform_n
)
{
	// ------------------------------ //
	// Sorting the scale coefficients //
	// ------------------------------ //
	
	void* d_temp_storage = NULL;
	size_t temp_storage  = 0;

	// this launch only decides how much temp_storage is needed for allocation to d_temp_storage
	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_assembled_solution.z0, 
		d_buffer_assembled_solution.z0, 
		d_assembled_solution.length
	) );

	d_temp_storage = malloc_device(temp_storage);

	if (non_uniform_n) {
		CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
		(
			d_temp_storage,
			temp_storage,
			d_morton_codes,
			d_sorted_morton_codes,
			d_assembled_solution.n0,
			d_buffer_assembled_solution.n0,
			d_assembled_solution.length
		));
	}


	// sorting the Morton codes is equivalent to ordering the scale coefficients according to
	// a z-order curve, please see: https://en.wikipedia.org/wiki/Z-order_curve
	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_assembled_solution.h, 
		d_buffer_assembled_solution.h, 
		d_assembled_solution.length
	) );

	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_assembled_solution.v, 
		d_buffer_assembled_solution.v, 
		d_assembled_solution.length
	) );

	//CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_morton_codes,
	//	d_sorted_morton_codes, 
	//	d_assembled_solution.q2, 
	//	d_buffer_assembled_solution.q2, 
	//	d_assembled_solution.length
	//) );

	//CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_morton_codes,
	//	d_sorted_morton_codes,
	//	d_assembled_solution.q3,
	//	d_buffer_assembled_solution.q3,
	//	d_assembled_solution.length
	//));

	//CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_morton_codes,
	//	d_sorted_morton_codes,
	//	d_assembled_solution.q4,
	//	d_buffer_assembled_solution.q4,
	//	d_assembled_solution.length
	//));

	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_assembled_solution.z0, 
		d_buffer_assembled_solution.z0, 
		d_assembled_solution.length
	) );

	CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes,
		d_assembled_solution.z1x,
		d_buffer_assembled_solution.z1x,
		d_assembled_solution.length
	));

	CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes,
		d_assembled_solution.z1y,
		d_buffer_assembled_solution.z1y,
		d_assembled_solution.length
	));
	
	//CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_morton_codes,
	//	d_sorted_morton_codes,
	//	d_assembled_solution.zxy,
	//	d_buffer_assembled_solution.zxy,
	//	d_assembled_solution.length
	//));

	free_device(d_temp_storage);

	// ------------------------------ //

	// ---------------------------------------------- //
	// Getting array with which to reverse z-ordering //
	// ---------------------------------------------- //

	d_temp_storage = NULL;
	temp_storage   = 0;

	// this launch only decides how much temp_storage is needed for allocation to d_temp_storage
	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_indices, 
		d_reverse_z_order,
		d_assembled_solution.length
	) );

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR( cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_morton_codes,
		d_sorted_morton_codes, 
		d_indices,
		d_reverse_z_order,
		d_assembled_solution.length
	) );

	free_device(d_temp_storage);

	// ---------------------------------------------- //
}