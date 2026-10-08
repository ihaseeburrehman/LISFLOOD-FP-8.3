#include "get_reg_tree.cuh"

__host__ void lis::cuda::acc_nugrid::get_reg_tree
(
	bool*         d_sig_details,
	int           lev
)
{
	for (int level = lev - 1; level > LVL_SINGLE_BLOCK; level--)
	{		
		index_1D prev_lvl_idx = get_lvl_idx(level - 1);
		index_1D curr_lvl_idx = get_lvl_idx(level);
		index_1D next_lvl_idx = get_lvl_idx(level + 1);
	    int      num_threads  = 1 << (2 * level);
		
		int num_blocks = get_num_blocks(num_threads, THREADS_PER_BLOCK_MRA);

		reg_new<false><<<num_blocks, THREADS_PER_BLOCK_MRA>>>
		(
			d_sig_details,
			level,
			prev_lvl_idx,
			curr_lvl_idx,
			next_lvl_idx,
			num_threads
		);
	}

	index_1D prev_lvl_idx = get_lvl_idx(LVL_SINGLE_BLOCK - 1);
	index_1D curr_lvl_idx = get_lvl_idx(LVL_SINGLE_BLOCK);
	index_1D next_lvl_idx = get_lvl_idx(LVL_SINGLE_BLOCK + 1);
    int      num_threads  = 1 << (2 * LVL_SINGLE_BLOCK);

	reg_new<true><<<1, THREADS_PER_BLOCK_MRA>>>
	(
		d_sig_details,
		LVL_SINGLE_BLOCK,
		prev_lvl_idx,
		curr_lvl_idx,
		next_lvl_idx,
		num_threads
	);
}