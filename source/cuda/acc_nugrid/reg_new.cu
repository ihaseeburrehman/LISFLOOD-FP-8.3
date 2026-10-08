#include "reg_new.cuh"

template <bool SINGLE_BLOCK>
__global__ void lis::cuda::acc_nugrid::reg_new
(
	bool*    d_sig_details,
	int      level,
	index_1D prev_lvl_idx,
	index_1D curr_lvl_idx,
	index_1D next_lvl_idx,
	int      num_threads
)
{
	__shared__ bool shared_sig_details[THREADS_PER_BLOCK_MRA];
	
	DetailChildren child_details;

	index_1D t_idx = threadIdx.x;
	index_1D idx   = blockIdx.x * blockDim.x + t_idx;

	if (idx >= num_threads) return;

	if (SINGLE_BLOCK)
	{
		index_1D parent_idx;
		index_1D child_idx = curr_lvl_idx + t_idx;
		
		shared_sig_details[t_idx] = d_sig_details[child_idx];

		__syncthreads();
		
		for (int lvl = LVL_SINGLE_BLOCK - 1; lvl >= 0; lvl--)
		{
			index_1D curr_lvl_idx_block = get_lvl_idx(lvl);
			int      num_threads        = 1 << (2 * lvl);

			parent_idx = curr_lvl_idx_block + t_idx;

			if (t_idx < num_threads)
			{
				child_details = get_child_details
				(
					shared_sig_details,
					4 * t_idx
				);
			}

			__syncthreads();

			if (t_idx < num_threads)
			{
				if (child_details.has_sig_detail()) d_sig_details[parent_idx] = SIGNIFICANT;

				shared_sig_details[t_idx] = child_details.has_sig_detail();
			}

			__syncthreads();
		}
	}
	else
	{
		index_1D g_idx = curr_lvl_idx + idx;

		shared_sig_details[t_idx] = d_sig_details[g_idx];

		__syncthreads();
		
		if ( t_idx >= (THREADS_PER_BLOCK_MRA / 4) ) return;

		index_1D t_idx_shifted = 4 * t_idx;
				 g_idx         = prev_lvl_idx + t_idx + blockIdx.x * (THREADS_PER_BLOCK_MRA / 4);

		child_details = get_child_details
		(
			shared_sig_details,
			t_idx_shifted
		);

		if ( child_details.has_sig_detail() ) d_sig_details[g_idx] = SIGNIFICANT;
	}	
}

inline void dummy_template_instantiator
(
	bool*    d_sig_details,
	int      level,
	index_1D prev_lvl_idx,
	index_1D curr_lvl_idx,
	index_1D next_lvl_idx,
	int      num_threads)
{
	lis::cuda::acc_nugrid::reg_new<true><<<1, 1>>>
	(
		d_sig_details, 
		level, 
		prev_lvl_idx,
		curr_lvl_idx, 
		next_lvl_idx, 
		num_threads
	);

	lis::cuda::acc_nugrid::reg_new<false><<<1, 1>>>
	(
		d_sig_details,
		level,
		prev_lvl_idx,
		curr_lvl_idx,
		next_lvl_idx,
		num_threads
	);
}