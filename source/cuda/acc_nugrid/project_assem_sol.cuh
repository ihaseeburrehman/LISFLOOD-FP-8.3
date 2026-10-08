#pragma once

#include "traverse_tree_of_sig_details.cuh"
#include "reverse_z_order_assembled_solution.cuh"
#include "write_reals_to_file.cuh"
#include "get_num_blocks.h"
#include "reinsert_assem_sol.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

__host__ AssembledSolution project_assem_sol
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
);

}
}
}