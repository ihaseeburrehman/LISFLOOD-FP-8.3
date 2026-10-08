#include "compaction.cuh"

void lis::cuda::acc_nugrid::compaction
(
	AssembledSolution& d_sol, 
	AssembledSolution& d_buf_sol, 
	CompactionFlags&   d_compaction_flags,
	int                num_finest_elems,
	bool               non_uniform_n,
	int                startfile
)
{
	void*  d_temp_storage = NULL;
	size_t temp_storage_bytes = 0;

	int* d_sol_len = (int*)malloc_device( sizeof(int) );
	int* h_sol_len = new int;

	CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	(
		d_temp_storage, 
		temp_storage_bytes, 
		d_buf_sol.z0, 
		d_compaction_flags.north_east, 
		d_sol.z0, 
		d_sol_len, 
		num_finest_elems
	) );

	d_temp_storage = malloc_device(temp_storage_bytes);

	//CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	//(
	//	d_temp_storage, 
	//	temp_storage_bytes, 
	//	d_buf_sol.h, 
	//	d_compaction_flags.north_east, 
	//	d_sol.h, 
	//	d_sol_len, 
	//	num_finest_elems
	//) );

	//CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	//(
	//	d_temp_storage, 
	//	temp_storage_bytes, 
	//	d_buf_sol.qx, 
	//	d_compaction_flags.north_east, 
	//	d_sol.qx, 
	//	d_sol_len, 
	//	num_finest_elems
	//) );

	//CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	//(
	//	d_temp_storage, 
	//	temp_storage_bytes, 
	//	d_buf_sol.qy, 
	//	d_compaction_flags.north_east, 
	//	d_sol.qy, 
	//	d_sol_len, 
	//	num_finest_elems
	//) );

	CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	(
		d_temp_storage, 
		temp_storage_bytes, 
		d_buf_sol.z0, 
		d_compaction_flags.north_east, 
		d_sol.z0, 
		d_sol_len, 
		num_finest_elems
	) );

	CHECK_CUDA_ERROR(cub::DeviceSelect::Flagged
	(
		d_temp_storage,
		temp_storage_bytes,
		d_buf_sol.z1x,
		d_compaction_flags.north_east,
		d_sol.z1x,
		d_sol_len,
		num_finest_elems
	));

	CHECK_CUDA_ERROR(cub::DeviceSelect::Flagged
	(
		d_temp_storage,
		temp_storage_bytes,
		d_buf_sol.z1y,
		d_compaction_flags.north_east,
		d_sol.z1y,
		d_sol_len,
		num_finest_elems
	));

	if (non_uniform_n) {
		CHECK_CUDA_ERROR(cub::DeviceSelect::Flagged
		(
			d_temp_storage,
			temp_storage_bytes,
			d_buf_sol.n0,
			d_compaction_flags.north_east,
			d_sol.n0,
			d_sol_len,
			num_finest_elems
		));
	}

	if (startfile) {
		CHECK_CUDA_ERROR(cub::DeviceSelect::Flagged
		(
			d_temp_storage,
			temp_storage_bytes,
			d_buf_sol.h,
			d_compaction_flags.north_east,
			d_sol.h,
			d_sol_len,
			num_finest_elems
		));
	}

	CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	(
		d_temp_storage, 
		temp_storage_bytes, 
		d_buf_sol.levels, 
		d_compaction_flags.north_east, 
		d_sol.levels, 
		d_sol_len, 
		num_finest_elems
	) );

	CHECK_CUDA_ERROR( cub::DeviceSelect::Flagged
	(
		d_temp_storage, 
		temp_storage_bytes, 
		d_buf_sol.act_idcs, 
		d_compaction_flags.north_east, 
		d_sol.act_idcs, 
		d_sol_len, 
		num_finest_elems
	) );

	copy_cuda
	(
		h_sol_len, 
		d_sol_len, 
		sizeof(int)
	);
	
	d_sol.length = *h_sol_len;

	free_device(d_temp_storage);
	free_device(d_sol_len);
	delete      h_sol_len;
}