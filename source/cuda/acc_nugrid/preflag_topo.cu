#include "preflag_topo.cuh"

__host__ void lis::cuda::acc_nugrid::preflag_topo
(
	ScaleCoefficients& d_scale_coeffs,
	Details&           d_details,
	bool*              d_sig_details,
	Maxes&             maxes, 
	NUMERIC_TYPE&              eps,
	int&               lev,
	bool               non_uniform_n,
	int                startfile
)
{
	for (int level = lev - 1; level >= 0; level--)
	{
		NUMERIC_TYPE epsilon_local = eps / ( 1 << (lev - level) );
		
	    int num_threads  = 1 << (2 * level);
		
		int num_blocks = get_num_blocks(num_threads, THREADS_PER_BLOCK_MRA);

		encode_and_thresh_topo<<<num_blocks, THREADS_PER_BLOCK_MRA>>>
		(
			d_scale_coeffs,
			d_details,
			d_sig_details, /// not used
			maxes,
			epsilon_local,
			level,
			non_uniform_n,
			startfile
		);
	}
}