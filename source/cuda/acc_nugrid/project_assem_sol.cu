#include "project_assem_sol.cuh"

__host__ lis::cuda::acc_nugrid::AssembledSolution lis::cuda::acc_nugrid::project_assem_sol
(
	const int&                  mesh_dim,
	bool*&                      d_sig_details,
	const ScaleCoefficients&    d_scale_coeffs,
	AssembledSolution           d_buf_assem_sol,
	const int&     lev,
	MortonCode*                 d_reverse_z_order,
	MortonCode*                 d_indices,
	AssembledSolution           d_assem_sol,
	AssembledSolution           d_plot_assem_sol
)
{
	int num_blocks_sol = get_num_blocks(d_assem_sol.length, THREADS_PER_BLOCK_MRA);

	reinsert_assem_sol<<<num_blocks_sol, THREADS_PER_BLOCK_MRA>>>
	(
		d_assem_sol, 
		d_scale_coeffs
	);

	int num_finest_elems      = mesh_dim * mesh_dim;
	int num_threads_traversal = num_finest_elems / 4;
	int num_blocks_traversal  = get_num_blocks(num_threads_traversal, THREADS_PER_BLOCK_MRA);
	
	traverse_tree_of_sig_details<<<num_blocks_traversal, THREADS_PER_BLOCK_MRA>>>
	(
		d_sig_details, 
		d_scale_coeffs,
		d_buf_assem_sol,
		num_threads_traversal,
		lev
	);

	reverse_z_order_assembled_solution
	(
		d_reverse_z_order,
		d_indices,
		d_buf_assem_sol,
		d_plot_assem_sol,
		num_finest_elems
	);

	return d_plot_assem_sol;
}