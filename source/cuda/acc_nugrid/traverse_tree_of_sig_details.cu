#include "traverse_tree_of_sig_details.cuh"

__global__ void lis::cuda::acc_nugrid::traverse_tree_of_sig_details
(
	bool*             d_sig_details,
	ScaleCoefficients d_scale_coeffs,
	AssembledSolution d_buf_assem_sol,
	int               num_threads,
	int  lev
)
{	
	__shared__ union
	{
		NUMERIC_TYPE     coeffs [4 * THREADS_PER_BLOCK_MRA];
		int      levels [4 * THREADS_PER_BLOCK_MRA];
		index_1D indices[4 * THREADS_PER_BLOCK_MRA];

	} shared;
	
	NUMERIC_TYPE     h_arr[4];
	NUMERIC_TYPE     v_arr[4];
	//NUMERIC_TYPE     qx_arr[4];
	//NUMERIC_TYPE     qy_arr[4];
	//NUMERIC_TYPE     z_arr[4];
	NUMERIC_TYPE     z0_arr[4];
	NUMERIC_TYPE     z1x_arr[4];
	NUMERIC_TYPE     z1y_arr[4];
//	NUMERIC_TYPE     zxy_arr[4];
	index_1D indices[4];
	int      levels[4];

	index_1D t_idx = threadIdx.x;
	index_1D idx   = blockIdx.x * blockDim.x + t_idx;

	int block_store_step = 3 * blockIdx.x * THREADS_PER_BLOCK_MRA;

	if (idx >= num_threads) return;

	index_1D  stack[52];
	index_1D *stackPtr = stack;
	*stackPtr++ = NULL;

	index_1D g_idx = 0;

	MortonCode fine_code = 4 * idx;

	int level = 0;

	do
	{
		index_1D curr_lvl_idx = get_lvl_idx(level);

		index_1D local_idx = g_idx - curr_lvl_idx;
		
		MortonCode current_code = local_idx;

		bool is_child = ( ( fine_code >> ( 2 * (lev - level) ) ) == current_code);

		if (is_child)
		{
			bool is_sig = d_sig_details[g_idx];

			if (!is_sig)
			{				
				NUMERIC_TYPE h = d_scale_coeffs.h[g_idx];
				NUMERIC_TYPE v = d_scale_coeffs.v[g_idx];
				//NUMERIC_TYPE qx  = d_scale_coeffs.qx[g_idx];
				//NUMERIC_TYPE qy  = d_scale_coeffs.qy[g_idx];
				//NUMERIC_TYPE z   = d_scale_coeffs.z[g_idx];

				NUMERIC_TYPE z0 = d_scale_coeffs.z0[g_idx];
				NUMERIC_TYPE z1x = d_scale_coeffs.z1x[g_idx];
				NUMERIC_TYPE z1y = d_scale_coeffs.z1y[g_idx];
//				NUMERIC_TYPE zxy = d_scale_coeffs.zxy[g_idx];

				#pragma unroll
				for (int i = 0; i < 4; i++)
				{
					h_arr[i] = h; // eta - z;
					v_arr[i] = v;
					//qx_arr[i]  = qx;
					//qy_arr[i]  = qy;
					//z_arr[i]   = z;
					z0_arr[i] = z0;
					z1x_arr[i] = z1x;
					z1y_arr[i] = z1y;
//					zxy_arr[i] = zxy;

					indices[i] = g_idx;
					levels[i]  = level;
				}

				goto store;
			}
			else
			{
				bool penultimate_level = (++level == lev);
				
				index_1D next_lvl_idx = get_lvl_idx(level);
				
				index_1D child_idx = next_lvl_idx + 4 * local_idx;

				if (!penultimate_level)
				{
					// get child indices and make index child_0 of current sub-element
					g_idx       = child_idx + 0;
					*stackPtr++ = child_idx + 1;
					*stackPtr++ = child_idx + 2;
					*stackPtr++ = child_idx + 3;
				}
				else
				{
					// reached penultimate level, add information to last level and exit
					
					#pragma unroll
					for (int i = 0; i < 4; i++)
					{
						// disrupted ordering of coefficients to ensure cache hit for z, as z is used twice
						// once to calculate h, and once when simply storing z
						//NUMERIC_TYPE eta = d_scale_coeffs.eta[child_idx + i];
						//NUMERIC_TYPE z   = d_scale_coeffs.z[child_idx + i];

						h_arr[i] = d_scale_coeffs.h[child_idx + i]; // eta - z;
						v_arr[i] = d_scale_coeffs.v[child_idx + i];
						//z_arr[i]   = z;
						//qx_arr[i]  = d_scale_coeffs.qx[child_idx + i];
						//qy_arr[i]  = d_scale_coeffs.qy[child_idx + i];
						z0_arr[i] = d_scale_coeffs.z0[child_idx + i];
						z1x_arr[i] = d_scale_coeffs.z1x[child_idx + i];
						z1y_arr[i] = d_scale_coeffs.z1y[child_idx + i];
//						zxy_arr[i] = d_scale_coeffs.zxy[child_idx + i];

						indices[i] = child_idx + i;
						levels[i]  = level;
					}
					
					goto store;
				}
			}
		}
		else
		{
			g_idx = *--stackPtr;
		}
	}
	while (NULL != g_idx);

	store:
	{
		// storing h
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = h_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.h[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing v
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = v_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.v[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		//// storing qy
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = qy_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.qy[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();

		// storing z0
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z0_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z0[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing z1x
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z1x_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z1x[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing z1y
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z1y_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z1y[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing zxy
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = zxy_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.zxy[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();
		
		// storing active indices
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.indices[4 * t_idx + i] = indices[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.act_idcs[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.indices[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();
		
		// storing levels
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.levels[4 * t_idx + i] = levels[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.levels[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.levels[t_idx + i * THREADS_PER_BLOCK_MRA];
	}
}

__global__ void lis::cuda::acc_nugrid::traverse_tree_of_sig_details_with_n
(
	bool* d_sig_details,
	ScaleCoefficients d_scale_coeffs,
	AssembledSolution d_buf_assem_sol,
	int               num_threads,
	int  lev
)
{
	__shared__ union
	{
		NUMERIC_TYPE     coeffs[4 * THREADS_PER_BLOCK_MRA];
		int      levels[4 * THREADS_PER_BLOCK_MRA];
		index_1D indices[4 * THREADS_PER_BLOCK_MRA];

	} shared;

	NUMERIC_TYPE     h_arr[4];
	NUMERIC_TYPE     v_arr[4];
	//NUMERIC_TYPE     qx_arr[4];
	//NUMERIC_TYPE     qy_arr[4];
	//NUMERIC_TYPE     z_arr[4];
	NUMERIC_TYPE     z0_arr[4];
	NUMERIC_TYPE     z1x_arr[4];
	NUMERIC_TYPE     z1y_arr[4];
	NUMERIC_TYPE     n0_arr[4];
	//	NUMERIC_TYPE     n1x_arr[4];
	//	NUMERIC_TYPE     n1y_arr[4];
	//	NUMERIC_TYPE     zxy_arr[4];
	index_1D indices[4];
	int      levels[4];

	index_1D t_idx = threadIdx.x;
	index_1D idx = blockIdx.x * blockDim.x + t_idx;

	int block_store_step = 3 * blockIdx.x * THREADS_PER_BLOCK_MRA;

	if (idx >= num_threads) return;

	index_1D  stack[52];
	index_1D* stackPtr = stack;
	*stackPtr++ = NULL;

	index_1D g_idx = 0;

	MortonCode fine_code = 4 * idx;

	int level = 0;

	do
	{
		index_1D curr_lvl_idx = get_lvl_idx(level);

		index_1D local_idx = g_idx - curr_lvl_idx;

		MortonCode current_code = local_idx;

		bool is_child = ((fine_code >> (2 * (lev - level))) == current_code);

		if (is_child)
		{
			bool is_sig = d_sig_details[g_idx];

			if (!is_sig)
			{
				NUMERIC_TYPE h = d_scale_coeffs.h[g_idx];
				NUMERIC_TYPE v = d_scale_coeffs.v[g_idx];
				//NUMERIC_TYPE qx  = d_scale_coeffs.qx[g_idx];
				//NUMERIC_TYPE qy  = d_scale_coeffs.qy[g_idx];
				//NUMERIC_TYPE z   = d_scale_coeffs.z[g_idx];

				NUMERIC_TYPE z0 = d_scale_coeffs.z0[g_idx];
				NUMERIC_TYPE z1x = d_scale_coeffs.z1x[g_idx];
				NUMERIC_TYPE z1y = d_scale_coeffs.z1y[g_idx];
				NUMERIC_TYPE n0 = d_scale_coeffs.n0[g_idx];
				//				NUMERIC_TYPE zxy = d_scale_coeffs.zxy[g_idx];

#pragma unroll
				for (int i = 0; i < 4; i++)
				{
					h_arr[i] = h; // eta - z;
					v_arr[i] = v; // eta - z;
					//qx_arr[i]  = qx;
					//qy_arr[i]  = qy;
					//z_arr[i]   = z;
					z0_arr[i] = z0;
					z1x_arr[i] = z1x;
					z1y_arr[i] = z1y;
					n0_arr[i] = n0;
					//					n1x_arr[i] = C(0.0);
					//					n1y_arr[i] = C(0.0);
					//					zxy_arr[i] = zxy;

					indices[i] = g_idx;
					levels[i] = level;
				}

				goto store;
			}
			else
			{
				bool penultimate_level = (++level == lev);

				index_1D next_lvl_idx = get_lvl_idx(level);

				index_1D child_idx = next_lvl_idx + 4 * local_idx;

				if (!penultimate_level)
				{
					// get child indices and make index child_0 of current sub-element
					g_idx = child_idx + 0;
					*stackPtr++ = child_idx + 1;
					*stackPtr++ = child_idx + 2;
					*stackPtr++ = child_idx + 3;
				}
				else
				{
					// reached penultimate level, add information to last level and exit

#pragma unroll
					for (int i = 0; i < 4; i++)
					{
						// disrupted ordering of coefficients to ensure cache hit for z, as z is used twice
						// once to calculate h, and once when simply storing z
						//NUMERIC_TYPE eta = d_scale_coeffs.eta[child_idx + i];
						//NUMERIC_TYPE z   = d_scale_coeffs.z[child_idx + i];

						h_arr[i] = d_scale_coeffs.h[child_idx + i]; // eta - z;
						v_arr[i] = d_scale_coeffs.v[child_idx + i]; // eta - z;
						//z_arr[i]   = z;
						//qx_arr[i]  = d_scale_coeffs.qx[child_idx + i];
						//qy_arr[i]  = d_scale_coeffs.qy[child_idx + i];
						z0_arr[i] = d_scale_coeffs.z0[child_idx + i];
						z1x_arr[i] = d_scale_coeffs.z1x[child_idx + i];
						z1y_arr[i] = d_scale_coeffs.z1y[child_idx + i];
						//						zxy_arr[i] = d_scale_coeffs.zxy[child_idx + i];
						n0_arr[i] = d_scale_coeffs.n0[child_idx + i];
						//						n1x_arr[i] = C(0.0);
						//						n1y_arr[i] = C(0.0);

						indices[i] = child_idx + i;
						levels[i] = level;
					}

					goto store;
				}
			}
		}
		else
		{
			g_idx = *--stackPtr;
		}
	} while (NULL != g_idx);

store:
	{
		// storing h
#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = h_arr[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.h[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing qx
		#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = v_arr[i];
		__syncthreads();
		#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.v[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		//// storing qy
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = qy_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.qy[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();

		// storing z0
#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z0_arr[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z0[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing z1x
#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z1x_arr[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z1x[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing z1y
#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = z1y_arr[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.z1y[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();




		// storing z0
#pragma unroll
		for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = n0_arr[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.n0[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		//// storing z1x
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = n1x_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.z1x[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();

		//// storing z1y
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = n1y_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.z1y[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();






		// storing zxy
		//#pragma unroll
		//for (int i = 0; i < 4; i++) shared.coeffs[4 * t_idx + i] = zxy_arr[i];
		//__syncthreads();
		//#pragma unroll
		//for (int i = 0; i < 4; i++) d_buf_assem_sol.zxy[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.coeffs[t_idx + i * THREADS_PER_BLOCK_MRA];
		//__syncthreads();

		// storing active indices
#pragma unroll
		for (int i = 0; i < 4; i++) shared.indices[4 * t_idx + i] = indices[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.act_idcs[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.indices[t_idx + i * THREADS_PER_BLOCK_MRA];
		__syncthreads();

		// storing levels
#pragma unroll
		for (int i = 0; i < 4; i++) shared.levels[4 * t_idx + i] = levels[i];
		__syncthreads();
#pragma unroll
		for (int i = 0; i < 4; i++) d_buf_assem_sol.levels[idx + i * THREADS_PER_BLOCK_MRA + block_store_step] = shared.levels[t_idx + i * THREADS_PER_BLOCK_MRA];
	}
}