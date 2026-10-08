#pragma once

#include "cuda_utils.cuh"
#include "../../lisflood.h"
#include "index_1D.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct NonUniformNeighbours
{
	NUMERIC_TYPE*     h_owner;
	NUMERIC_TYPE*     h_nghbr;
	NUMERIC_TYPE*     z_owner;
	NUMERIC_TYPE*     z_nghbr;
	NUMERIC_TYPE*     n_owner;
	NUMERIC_TYPE*     n_nghbr;
	NUMERIC_TYPE*     q;
	NUMERIC_TYPE*     z;
	NUMERIC_TYPE*     v;
	int*      nghbr_elem_idcs;
	index_1D* nghbr_act_idcs;
	int*      nghbr_flags;
	int*      owner_elem_idcs;
	index_1D* owner_act_idcs;
	int*      owner_counts;
	int*      owner_flags;
	NUMERIC_TYPE*     dx;
	NUMERIC_TYPE*     dx_sum;
	int*      dirs;
	int*      itrface_counts;
	int*      cumu_itrface_counts;
	int*      sum_flags;
	int       length;
	bool      is_copy;

	NonUniformNeighbours(const int& num_nghbrs, bool non_uniform_n)
	{
		size_t bytes_real = num_nghbrs * sizeof(NUMERIC_TYPE);
		size_t bytes_int  = num_nghbrs * sizeof(index_1D); // MKS: or int???

		h_owner             = (NUMERIC_TYPE*)malloc_device(bytes_real);
		h_nghbr             = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z_owner             = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z_nghbr             = (NUMERIC_TYPE*)malloc_device(bytes_real);
		q                   = (NUMERIC_TYPE*)malloc_device(bytes_real);
		z                   = (NUMERIC_TYPE*)malloc_device(bytes_real);
		v                   = (NUMERIC_TYPE*)malloc_device(bytes_real);

		n_owner = (non_uniform_n) ? (NUMERIC_TYPE*)malloc_device(bytes_real) : nullptr;
		n_nghbr = (non_uniform_n) ? (NUMERIC_TYPE*)malloc_device(bytes_real) : nullptr;

		nghbr_elem_idcs     = (int*)malloc_device(bytes_int);
		nghbr_act_idcs      = (index_1D*)malloc_device(bytes_int);
		nghbr_flags         = (int*)malloc_device(bytes_int);
		owner_elem_idcs     = (int*)malloc_device(bytes_int);
		owner_act_idcs      = (index_1D*)malloc_device(bytes_int);
		owner_counts        = (int*)malloc_device(bytes_int);
		owner_flags         = (int*)malloc_device(bytes_int);
		dx                  = (NUMERIC_TYPE*)malloc_device(bytes_real);
		dx_sum              = (NUMERIC_TYPE*)malloc_device(bytes_real);
		dirs                = (int*)malloc_device(bytes_int);
		itrface_counts      = (int*)malloc_device(bytes_int);
		cumu_itrface_counts = (int*)malloc_device(bytes_int);
		sum_flags           = (int*)malloc_device(bytes_int);
		length              = num_nghbrs;
		is_copy             = false;
	}

	NonUniformNeighbours(const NonUniformNeighbours& original) { *this = original; is_copy = true; }

	~NonUniformNeighbours()
	{
		if (!is_copy)
		{
			free_device(h_owner);
			free_device(h_nghbr);
			free_device(z_owner);
			free_device(z_nghbr);
			free_device(n_owner);
			free_device(n_nghbr);
			free_device(q);
			free_device(z);
			free_device(v);
			free_device(nghbr_elem_idcs);
			free_device(nghbr_act_idcs);
			free_device(nghbr_flags);
			free_device(owner_elem_idcs);
			free_device(owner_act_idcs);
			free_device(dx);
			free_device(dx_sum);
			free_device(dirs);
			free_device(owner_counts);
			free_device(owner_flags);
			free_device(itrface_counts);
			free_device(cumu_itrface_counts);
			free_device(sum_flags);
		}
	}

} NonUniformNeighbours;

}
}
}