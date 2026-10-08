#include "init_interfaces.cuh"
#include <helper_cuda.h> // MKS
#include <cub/cub.cuh> // MKS

__host__ lis::cuda::acc_nugrid::NonUniformInterfaces lis::cuda::acc_nugrid::init_interfaces
(
	NonUniformNeighbours d_non_uniform_nghbrs
)
{
	NUMERIC_TYPE* h_num_interfaces = new NUMERIC_TYPE;
	NUMERIC_TYPE* d_num_interfaces = (NUMERIC_TYPE*)malloc_device(sizeof(NUMERIC_TYPE));

	void* d_temp_storage = NULL;
	size_t temp_storage = 0;

	CHECK_CUDA_ERROR(cub::DeviceReduce::Sum
	(
		d_temp_storage,
		temp_storage,
		d_non_uniform_nghbrs.itrface_counts,
		d_num_interfaces,
		d_non_uniform_nghbrs.length
	));

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR(cub::DeviceReduce::Sum
	(
		d_temp_storage,
		temp_storage,
		d_non_uniform_nghbrs.itrface_counts,
		d_num_interfaces,
		d_non_uniform_nghbrs.length
	));

	copy_cuda
	(
		h_num_interfaces,
		d_num_interfaces,
		sizeof(NUMERIC_TYPE)
	);

	int num_interfaces = *h_num_interfaces;

	free_device(d_num_interfaces);
	free_device(d_temp_storage);
	delete      h_num_interfaces;

	// FINDING CUMULATIVE INTERFACE COUNTS //

	d_temp_storage = NULL;
	temp_storage = 0;

	CHECK_CUDA_ERROR(cub::DeviceScan::ExclusiveSum
	(
		d_temp_storage,
		temp_storage,
		d_non_uniform_nghbrs.itrface_counts,
		d_non_uniform_nghbrs.cumu_itrface_counts,
		d_non_uniform_nghbrs.length
	));

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR(cub::DeviceScan::ExclusiveSum
	(
		d_temp_storage,
		temp_storage,
		d_non_uniform_nghbrs.itrface_counts,
		d_non_uniform_nghbrs.cumu_itrface_counts,
		d_non_uniform_nghbrs.length
	));

	free_device(d_temp_storage);

	// ------------------------------- //

	return NonUniformInterfaces(num_interfaces);
}