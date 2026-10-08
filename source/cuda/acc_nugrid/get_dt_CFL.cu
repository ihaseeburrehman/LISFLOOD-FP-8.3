#include "get_dt_CFL.cuh"

__host__ NUMERIC_TYPE lis::cuda::acc_nugrid::get_dt_CFL
(
	NUMERIC_TYPE*&     d_dt_CFL,
	const int& sol_len
)
{
	NUMERIC_TYPE* h_min_out = new NUMERIC_TYPE;
	NUMERIC_TYPE* d_min_out = (NUMERIC_TYPE*)malloc_device(sizeof(NUMERIC_TYPE));

	void*  d_temp_storage = NULL;
	size_t temp_storage  = 0;

	CHECK_CUDA_ERROR( cub::DeviceReduce::Min
	(
		d_temp_storage,
		temp_storage,
		d_dt_CFL,
		d_min_out,
		sol_len
	) );

	d_temp_storage = malloc_device(temp_storage);

	CHECK_CUDA_ERROR( cub::DeviceReduce::Min
	(
		d_temp_storage,
		temp_storage,
		d_dt_CFL,
		d_min_out,
		sol_len
	) );

	copy_cuda
	(
		h_min_out, 
		d_min_out, 
		sizeof(NUMERIC_TYPE)
	);

	NUMERIC_TYPE dt_min = *h_min_out;

	free_device(d_min_out);
	free_device(d_temp_storage);
	delete h_min_out;

	return dt_min;
}