#include "sort_nghbr_assem_sol_row_major.cuh"
#include <helper_cuda.h> // MKS
#include <cub/cub.cuh> // MKS

__host__ void lis::cuda::acc_nugrid::sort_nghbr_assem_sol_row_major
(
	MortonCode*        d_reverse_z_order,
	MortonCode*        d_indices,
	AssembledSolution& d_buf_assem_sol,
	AssembledSolution& d_nghbr_assem_sol,
	int                array_length
)
{
	void* d_temp_storage = NULL;
	size_t temp_storage = 0;

	// launch to decides how much temp_storage is needed for allocation to d_temp_storage
	CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_reverse_z_order,
		d_indices,
		d_buf_assem_sol.act_idcs,
		d_nghbr_assem_sol.act_idcs,
		array_length
	));

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_reverse_z_order,
		d_indices,
		d_buf_assem_sol.act_idcs,
		d_nghbr_assem_sol.act_idcs,
		array_length
	));

	CHECK_CUDA_ERROR(cub::DeviceRadixSort::SortPairs
	(
		d_temp_storage,
		temp_storage,
		d_reverse_z_order,
		d_indices,
		d_buf_assem_sol.levels,
		d_nghbr_assem_sol.levels,
		array_length
	));

	free_device(d_temp_storage);
}