#pragma once
#include "stats.h"
#include "cuda_flow.cuh"
#include "../lisflood.h"

namespace lis
{
namespace cuda
{

typedef struct TimeSeries
{
	int count;
	NUMERIC_TYPE* __restrict__ time;
	NUMERIC_TYPE* __restrict__ value;
	NUMERIC_TYPE current_value;
} TimeSeries;

typedef struct BoundaryConditions
{
	// boundary conditions
	int BC_count;
	ESourceType* __restrict__ BC_type;
	NUMERIC_TYPE* __restrict__ BC_value;
	TimeSeries** __restrict__ BC_time_series;

	// point sources
	int PS_count;
	ESourceType* __restrict__ PS_type;
	NUMERIC_TYPE* __restrict__ PS_value;
	TimeSeries** __restrict__ PS_time_series;
	int* __restrict__ PS_idx;

	TimeSeries* __restrict__ all_time_series;
	int time_series_count;
} BoundaryConditions;

struct Boundary
{
	static void initialise
	(
		BoundaryConditions& d_dst,
		BoundCs& h_src,
		int pitch,
		int offset
	);

	static void update_time_series
	(
		NUMERIC_TYPE t
	);

	static void update_point_sources
	(
	 	NUMERIC_TYPE* H,
		NUMERIC_TYPE* DEM,
		MassStats* mass_stats
	);

	static void update_H
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* DEM,
		NUMERIC_TYPE* Qx,
		NUMERIC_TYPE* Qy,
		MassStats* mass_stats,
		dim3 grid_size,
		int drycheck
	);

	static void drain_nodata_water
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* DEM,
		MassStats* mass_stats,
		dim3 grid_size
	);

	static void drain_nodata_waterACC 
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* DEM,
		MassStats* mass_stats,
		dim3 grid_size
	);

	static void free_device
	(
	 	BoundaryConditions& boundaries
	);

	__device__ static FlowVector outside_x
	(
		FlowVector U_const,
		FlowVector U_inside,
		int bc_i,
		int HU_sign,
		NUMERIC_TYPE Zstar
	);

	__device__ static FlowVector outside_y
	(
		FlowVector U_const,
		FlowVector U_inside,
		int bc_i,
		int HU_sign,
		NUMERIC_TYPE Zstar
	);

	__device__ static FlowVector inside_x
	(
		FlowVector U_outside,
		FlowVector U_const,
		FlowVector U_inside,
		int bc_i
	);

	__device__ static FlowVector inside_y
	(
		FlowVector U_outside,
		FlowVector U_const,
		FlowVector U_inside,
		int bc_i
	);

	__device__ static int index_w
	(
		const int i,
		const int j
	);

	__device__ static int index_e
	(
		const int i,
		const int j
	);

	__device__ static int index_n
	(
		const int i,
		const int j
	);

	__device__ static int index_s
	(
		const int i,
		const int j
	);
	
	__device__ static int index_w_ACC 
	(
		const int i,
		const int j
	);

	__device__ static int index_e_ACC 
	(
		const int i,
		const int j
	);

	__device__ static int index_n_ACC 
	(
		const int i,
		const int j
	);

	__device__ static int index_s_ACC 
	(
		const int i,
		const int j
	);	

private:
	static void initialise_all_time_series
	(
		BoundaryConditions& dst,
		BoundCs& src
	);

	static void initialise_BC
	(
		BoundaryConditions& dst,
		BoundCs& src
	);

	static void initialise_PS
	(
		BoundaryConditions& dst,
		BoundCs& src,
		int pitch,
		int offset
	);
};

}
}
