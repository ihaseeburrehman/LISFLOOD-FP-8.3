#include "get_max_scale_coefficients.cuh"
#include <helper_cuda.h> // MKS
#include <cub/cub.cuh> // MKS

__host__ lis::cuda::acc_nugrid::Maxes lis::cuda::acc_nugrid::get_max_scale_coefficients
(
	AssembledSolution& d_assem_sol
)
{
	// Variables for maxes //
	Maxes maxes = { C(1.0) };

	NUMERIC_TYPE* h_max_out = new NUMERIC_TYPE;
	NUMERIC_TYPE* d_max_out = (NUMERIC_TYPE*)malloc_device( sizeof(NUMERIC_TYPE) );

	// --------------------//

	// Allocating memory to find maxes //

	void* d_temp_storage  = NULL;
	size_t  temp_storage  = 0;

	CHECK_CUDA_ERROR( cub::DeviceReduce::Max
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.z0, /// in
		d_max_out, /// out
		d_assem_sol.length
	) );

	d_temp_storage = malloc_device(temp_storage);

	// ------------------------------- //

	// Finding maxes //

	// eta
	//int num_blocks = get_num_blocks(d_assem_sol.length, THREADS_PER_BLOCK_MRA);
	//
	//init_eta_temp<<<num_blocks, THREADS_PER_BLOCK_MRA>>> /// adds h+z for the finest grid
	//(
	//	d_assem_sol, 
	//	d_eta_temp
	//);

	//CHECK_CUDA_ERROR( cub::DeviceReduce::Max
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_eta_temp,
	//	d_max_out, /// FMAX of eta on finest
	//	d_assem_sol.length
	//) );

	//copy_cuda
	//(
	//	h_max_out,
	//	d_max_out,
	//	sizeof(NUMERIC_TYPE)
	//);

	//maxes.eta = FMAX( *h_max_out, C(1.0) ); /// FMAX (eta_max, 1.0)

	// qx
	//CHECK_CUDA_ERROR( cub::DeviceReduce::Reduce
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_assem_sol.qx,
	//	d_max_out,
	//	d_assem_sol.length,
	//	abs_max,
	//	C(0.0)
	//) );

	//copy_cuda
	//(
	//	h_max_out,
	//	d_max_out,
	//	sizeof(NUMERIC_TYPE)
	//);

//	maxes.qx = FMAX( *h_max_out, C(1.0) );

	// qy
	//CHECK_CUDA_ERROR( cub::DeviceReduce::Reduce
	//(
	//	d_temp_storage,
	//	temp_storage,
	//	d_assem_sol.qy,
	//	d_max_out,
	//	d_assem_sol.length,
	//	abs_max,
	//	C(0.0)
	//) );

	//copy_cuda
	//(
	//	h_max_out,
	//	d_max_out,
	//	sizeof(NUMERIC_TYPE)
	//);

	//maxes.qy = FMAX( *h_max_out, C(1.0) );

	// z
	CHECK_CUDA_ERROR( cub::DeviceReduce::Max
	(
		d_temp_storage,
		temp_storage,
		d_assem_sol.z0,
		d_max_out,
		d_assem_sol.length
	) );

	copy_cuda
	(
		h_max_out,
		d_max_out,
		sizeof(NUMERIC_TYPE)
	);

	maxes.z = FMAX( *h_max_out, C(1.0) );

	// ------------- //

	free_device(d_max_out);
	free_device(d_temp_storage);
	delete h_max_out;

	return maxes;
}