#include "get_num_blocks.h"

int lis::cuda::acc_nugrid::get_num_blocks
(
	int num_threads, 
	int threads_per_block
)
{
	return num_threads / threads_per_block + (num_threads % threads_per_block != 0);
}