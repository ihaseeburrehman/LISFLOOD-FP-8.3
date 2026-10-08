#include "cuda_atomic.cuh"
#include "cuda_boundary.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"

namespace lis
{
namespace cuda
{

__device__ NUMERIC_TYPE linear_interpolate
(
	TimeSeries& time_series,
	NUMERIC_TYPE t
)
{
	NUMERIC_TYPE* times = time_series.time;
	NUMERIC_TYPE* values = time_series.value;

	if (t < times[0]) return values[0];

	for (int i=1; i < time_series.count; i++)
	{
		if (times[i - 1] <= t && times[i] > t)
		{
			NUMERIC_TYPE dt, a1, a2;
			dt = times[i] - times[i - 1];
			a1 = (t - times[i - 1]) / dt;
			a2 = C(1.0) - a1;

			return a1*values[i] + a2*values[i - 1];
		}
	}
	
	return values[time_series.count - 1];
}

__global__ void update_time_series
(
	NUMERIC_TYPE t
)
{
	for (int i=0; i < cuda::boundaries.time_series_count; i++)
	{
		TimeSeries& time_series = cuda::boundaries.all_time_series[i];
		time_series.current_value = linear_interpolate(time_series, t);
	}
}

__global__ void update_time_varying_boundary_conditions()
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;

	for (int i=global_i; i<cuda::boundaries.BC_count; i+=blockDim.y*gridDim.x)
	{
		if (cuda::boundaries.BC_type[i] == HVAR3
				|| cuda::boundaries.BC_type[i] == QVAR5)
		{
			cuda::boundaries.BC_value[i] =
				cuda::boundaries.BC_time_series[i]->current_value;
		}
	}
}

__global__ void update_time_varying_point_sources()
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;

	for (int i=global_i; i<cuda::boundaries.PS_count; i+=blockDim.y*gridDim.x)
	{
		if (cuda::boundaries.PS_type[i] == HVAR3
				|| cuda::boundaries.PS_type[i] == QVAR5)
		{
			cuda::boundaries.PS_value[i] =
				cuda::boundaries.PS_time_series[i]->current_value;
		}
	}
}

__device__ void accumulate_point_source_stats
(
	MassStats* mass_stats,
	NUMERIC_TYPE discharge
)
{
	if (discharge > C(0.0))
	{
		atomicAdd(&(mass_stats->in), discharge);
	}
	else
	{
		atomicAdd(&(mass_stats->out), -discharge);
	}
}

__global__ void update_point_sources
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;

	for (int i=global_i; i<cuda::boundaries.PS_count; i+=blockDim.x*gridDim.x)
	{
		NUMERIC_TYPE& Hvalue = H[cuda::boundaries.PS_idx[i]];

		switch (cuda::boundaries.PS_type[i])
		{
		case HFIX2:
		case HVAR3:
			{
				NUMERIC_TYPE Z = DEM[cuda::boundaries.PS_idx[i]];
				NUMERIC_TYPE H_new = FMIN(C(0.0),
						cuda::boundaries.PS_value[i] - Z);
				NUMERIC_TYPE discharge = (H_new - Hvalue)
					* cuda::geometry.dx*cuda::geometry.dy / cuda::dt;
				Hvalue = H_new;	

				accumulate_point_source_stats(mass_stats, discharge);
			}
			break;
		case QFIX4:
		case QVAR5:
			{
				NUMERIC_TYPE discharge = cuda::boundaries.PS_value[i];
				Hvalue += discharge * cuda::dt / cuda::geometry.dx;

				accumulate_point_source_stats(mass_stats,
						discharge * cuda::geometry.dx);
			}
			break;
		}
	}
}

__global__ void update_point_sources_Q
(
	NUMERIC_TYPE* H,
	MassStats* mass_stats
)
{
	int global_i = blockIdx.x * blockDim.x + threadIdx.x;

	for (int i = global_i; i < cuda::boundaries.PS_count; i += blockDim.x * gridDim.x)
	{
		NUMERIC_TYPE& Hvalue = H[cuda::boundaries.PS_idx[i]];

		switch (cuda::boundaries.PS_type[i])
		{
		case QFIX4:
		{
			NUMERIC_TYPE discharge = cuda::boundaries.PS_value[i]; //// like qtmp in UpdateH
			Hvalue += discharge * cuda::dt / cuda::geometry.dx;

			accumulate_point_source_stats(mass_stats, //// mass_stats->in is Qpoint_pos
				discharge * cuda::geometry.dx);   //// mass_stats->out is Qpoint_neg
		}
		break;
		case QVAR5:
		{
			NUMERIC_TYPE discharge = cuda::boundaries.PS_value[i]; //// like qtmp in UpdateH
			Hvalue += discharge * cuda::dt / cuda::geometry.dx;

			accumulate_point_source_stats(mass_stats, //// mass_stats->in is Qpoint_pos
				discharge * cuda::geometry.dx);   //// mass_stats->out is Qpoint_neg

		}
		break;
		}
	}
}

__global__ void update_point_sources_H
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats
)
{
	int global_i = blockIdx.x * blockDim.x + threadIdx.x;

	for (int i = global_i; i < cuda::boundaries.PS_count; i += blockDim.x * gridDim.x)
	{
		NUMERIC_TYPE& Hvalue = H[cuda::boundaries.PS_idx[i]];

		switch (cuda::boundaries.PS_type[i])
		{
		case HFIX2:
		{
			NUMERIC_TYPE H_new = cuda::boundaries.PS_value[i] - DEM[cuda::boundaries.PS_idx[i]]; //// H_new is himp of acc cpu

			if (H_new < C(0.0)) H_new = C(0.0);

			NUMERIC_TYPE discharge = (H_new - Hvalue) * cuda::geometry.dx * cuda::geometry.dx / cuda::dt;
			Hvalue = H_new;

			accumulate_point_source_stats(mass_stats, discharge);
		}
		break;
		case HVAR3:
		{
			NUMERIC_TYPE H_new = cuda::boundaries.PS_value[i] - DEM[cuda::boundaries.PS_idx[i]]; //// H_new is himp of acc cpu

			if (H_new < C(0.0)) H_new = C(0.0);

			NUMERIC_TYPE discharge = (H_new - Hvalue) * cuda::geometry.dx * cuda::geometry.dx / cuda::dt;
			Hvalue = H_new;

			accumulate_point_source_stats(mass_stats, discharge);
		}
		break;
		}
	}
}


__global__ void update_H
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* Qx,
	NUMERIC_TYPE* Qy,
	int drycheck
)
{
	NUMERIC_TYPE dV, cv, WDweight;

	int global_i = blockIdx.x * blockDim.x + threadIdx.x;
	int global_j = blockIdx.y * blockDim.y + threadIdx.y;

	for (int j = global_j; j < cuda::geometry.ysz; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::geometry.xsz; i += blockDim.x * gridDim.x)
		{

			NUMERIC_TYPE qxptr0 = Qx[i + j * (cuda::geometry.xsz + 1)];
			NUMERIC_TYPE qxptr1 = Qx[i + j * (cuda::geometry.xsz + 1) + 1];
			NUMERIC_TYPE qyptr0 = Qy[i + j * (cuda::geometry.xsz + 1)];
			NUMERIC_TYPE qyptr1 = Qy[i + (j + 1) * (cuda::geometry.xsz + 1)];
			NUMERIC_TYPE& hptr = H[i + j * (cuda::geometry.xsz)];
			
			dV = cuda::dt * (qxptr0 - qxptr1 + qyptr0 - qyptr1);

			cv = hptr * (cuda::geometry.dx * cuda::geometry.dx);

			if (drycheck) {

				if (cv + dV < C(0.0))
				{
					WDweight = -cv / dV; 
					if (qxptr0 < C(0.0)) qxptr0 *= WDweight;
					if (qxptr1 > C(0.0)) qxptr1 *= WDweight;
					if (qyptr0 < C(0.0)) qyptr0 *= WDweight;
					if (qyptr1 > C(0.0)) qyptr1 *= WDweight;
				}
			}

			dV = cuda::dt * (qxptr0 - qxptr1 + qyptr0 - qyptr1);
			hptr += dV / (cuda::geometry.dx * cuda::geometry.dx);
			if (hptr < C(0.0)) hptr = C(0.0);

		}
	}
}
__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
drain_nodata_water
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{

			NUMERIC_TYPE Z = DEM[j*cuda::pitch + i];
			if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6))
			{
				NUMERIC_TYPE& Hval = H[j*cuda::pitch + i];
				NUMERIC_TYPE Hold = Hval;
				Hval = C(0.0);

				NUMERIC_TYPE Q = -Hold * cuda::geometry.dx * cuda::geometry.dy
					/ cuda::dt;
				accumulate_point_source_stats(mass_stats, Q);
			}
		}
	}
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
drain_nodata_waterACC
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::geometry.xsz; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE Z = DEM[j*cuda::geometry.xsz + i];
			if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6))
			{
				NUMERIC_TYPE& Hval = H[j*cuda::pitch + i];
				NUMERIC_TYPE Hold = Hval;
				Hval = C(0.0);

				NUMERIC_TYPE Q = -Hold * cuda::geometry.dx * cuda::geometry.dy
					/ cuda::dt;
				accumulate_point_source_stats(mass_stats, Q);
			}
		}
	}
}			   
}
}

__device__ lis::cuda::FlowVector lis::cuda::Boundary::outside_x
(
	FlowVector U_const,
	FlowVector U_inside,
	int bc_i,
	int HU_sign,
	NUMERIC_TYPE Zstar
)
{
	FlowVector U_outside = { U_const.H, U_const.HU, U_const.HV }; 

	switch (cuda::boundaries.BC_type[bc_i])
	{
	case FREE1:
		break;
	case HFIX2:
	case HVAR3:
		U_outside.H = FMAX(C(0.0), cuda::boundaries.BC_value[bc_i] - Zstar);
		break;
	case QFIX4:
	case QVAR5:
		U_outside.HU = HU_sign * cuda::boundaries.BC_value[bc_i];
		U_outside.HV = C(0.0);
		if (FABS(U_outside.HU) > C(1e-6))
		{
			U_outside.H = FMAX(U_const.H,
					C(1.1)*cuda::solver_params.DepthThresh);
		}
		break;
	case NONE0:
	default:
		U_outside.HU = -U_const.HU;
		break;
	}

	return U_outside;
}

__device__ lis::cuda::FlowVector lis::cuda::Boundary::outside_y
(
	FlowVector U_const,
	FlowVector U_inside,
	int bc_i,
	int HU_sign,
	NUMERIC_TYPE Zstar
)
{
	FlowVector U_const_rotated = { U_const.H, U_const.HV, U_const.HU };
	FlowVector U_inside_rotated = { U_inside.H, U_inside.HV, U_inside.HU };
	FlowVector U_outside = Boundary::outside_x(U_const_rotated,
			U_inside_rotated, bc_i, HU_sign, Zstar);

	return { U_outside.H, U_outside.HV, U_outside.HU };
}

__device__ lis::cuda::FlowVector lis::cuda::Boundary::inside_x
(
	FlowVector U_outside,
	FlowVector U_const,
	FlowVector U_inside,
	int bc_i
)
{
	switch (cuda::boundaries.BC_type[bc_i])
	{
	case FREE1:
		return U_inside;
	case HFIX2:
	case HVAR3:
		U_const.H = U_outside.H;
		break;
	case QFIX4:
	case QVAR5:
		return { U_outside.H, U_outside.HU, U_const.HV };
	}

	return U_const;
}

__device__ lis::cuda::FlowVector lis::cuda::Boundary::inside_y
(
	FlowVector U_outside,
	FlowVector U_const,
	FlowVector U_inside,
	int bc_i
)
{
	FlowVector U_outside_rotated = { U_outside.H, U_outside.HV, U_outside.HU };
	FlowVector U_const_rotated = { U_const.H, U_const.HV, U_const.HU };
	FlowVector U_inside_rotated = { U_inside.H, U_inside.HV, U_inside.HU };
	U_inside = Boundary::inside_x(U_outside_rotated, U_const_rotated, U_inside_rotated, bc_i);
	return { U_inside.H, U_inside.HV, U_inside.HU };
}

__device__ int lis::cuda::Boundary::index_w
(
	const int i,
	const int j
)
{
	return 2*cuda::geometry.xsz + 2*cuda::geometry.ysz - j;
}

__device__ int lis::cuda::Boundary::index_e
(
	const int i,
	const int j
)
{
	return cuda::geometry.xsz + j-1;
}

__device__ int lis::cuda::Boundary::index_n
(
	const int i,
	const int j
)
{
	return i-1;
}

__device__ int lis::cuda::Boundary::index_s
(
	const int i,
	const int j
)
{
	return 2*cuda::geometry.xsz + cuda::geometry.ysz - i;
}

__device__ int lis::cuda::Boundary::index_w_ACC
(
	const int i,
	const int j
)
{
	return 2*cuda::geometry.xsz + 2*cuda::geometry.ysz - j - 1; 
}

__device__ int lis::cuda::Boundary::index_e_ACC
(
	const int i,
	const int j
)
{
	return cuda::geometry.xsz + j; 
}

__device__ int lis::cuda::Boundary::index_n_ACC
(
	const int i,
	const int j
)
{
	return i; 
}

__device__ int lis::cuda::Boundary::index_s_ACC
(
	const int i,
	const int j
)
{
	return 2*cuda::geometry.xsz + cuda::geometry.ysz - i - 1;
}											   
void lis::cuda::Boundary::initialise
(
	BoundaryConditions& d_dst,
	BoundCs& h_src,
	int pitch,
	int offset
)
{
	BoundaryConditions temp;

	initialise_all_time_series(temp, h_src);
	initialise_BC(temp, h_src);
	initialise_PS(temp, h_src, pitch, offset);
	
	cuda::copy_to_symbol(d_dst, &temp, sizeof(BoundaryConditions));
}

void lis::cuda::Boundary::initialise_all_time_series
(
	BoundaryConditions& dst,
	BoundCs& src
)
{
	dst.time_series_count = src.allTimeSeries.size();

	dst.all_time_series = static_cast<TimeSeries*>(cuda::malloc_device(
			dst.time_series_count*sizeof(TimeSeries)));

	TimeSeries* temp_all_time_series = new TimeSeries[dst.time_series_count];
	for (int i=0; i<dst.time_series_count; i++)
	{
		int count = src.allTimeSeries[i].count;
		temp_all_time_series[i].count = count;

		temp_all_time_series[i].time = static_cast<NUMERIC_TYPE*>(
				cuda::malloc_device(count*sizeof(NUMERIC_TYPE)));
		temp_all_time_series[i].value = static_cast<NUMERIC_TYPE*>(
				cuda::malloc_device(count*sizeof(NUMERIC_TYPE)));

		cuda::copy(temp_all_time_series[i].time, src.allTimeSeries[i].time,
				count*sizeof(NUMERIC_TYPE));
		cuda::copy(temp_all_time_series[i].value, src.allTimeSeries[i].value,
				count*sizeof(NUMERIC_TYPE));
	}
	cuda::copy(dst.all_time_series, temp_all_time_series,
			dst.time_series_count*sizeof(TimeSeries));
	delete[] temp_all_time_series;
}

void lis::cuda::Boundary::initialise_BC
(
	BoundaryConditions& dst,
	BoundCs& src
)
{
	dst.BC_count = src.numBCs;
	dst.BC_type = static_cast<ESourceType*>(
			cuda::malloc_device(src.numBCs*sizeof(ESourceType)));
	cuda::copy(dst.BC_type, src.BC_Ident, src.numBCs*sizeof(ESourceType));

	dst.BC_value = static_cast<NUMERIC_TYPE*>(
			cuda::malloc_device(src.numBCs*sizeof(NUMERIC_TYPE)));
	cuda::copy(dst.BC_value, src.BC_Val, src.numBCs*sizeof(NUMERIC_TYPE));

	dst.BC_time_series = static_cast<TimeSeries**>(
			cuda::malloc_device(src.numBCs*sizeof(TimeSeries*)));
	TimeSeries** time_series = new TimeSeries*[src.numBCs];

	for (int bc_i=0; bc_i<src.numBCs; bc_i++)
	{
		if (src.BC_Ident[bc_i] == HVAR3 || src.BC_Ident[bc_i] == QVAR5)
		{
			
			::TimeSeries* target = src.BC_TimeSeries[bc_i];

			for (std::vector<::TimeSeries>::size_type ts_i=0,
					e=src.allTimeSeries.size(); ts_i!=e; ts_i++)
			{
				::TimeSeries* candidate = &(src.allTimeSeries[ts_i]);
				if (candidate->value == target->value)
				{
					time_series[bc_i] = &(dst.all_time_series[ts_i]);
						//+ ts_i*sizeof(TimeSeries*);
				}
			}
		}
	}
	cuda::copy(dst.BC_time_series,
			time_series, src.numBCs*sizeof(TimeSeries*));
	delete[] time_series;
}

void lis::cuda::Boundary::initialise_PS
(
	BoundaryConditions& dst,
	BoundCs& src,
	int pitch,
	int offset
)
{
	int numPS = FMAX(0, src.numPS);
	dst.PS_count = numPS;

	dst.PS_type = static_cast<ESourceType*>(
			cuda::malloc_device(numPS*sizeof(ESourceType)));
	cuda::copy(dst.PS_type, src.PS_Ident, numPS*sizeof(ESourceType));

	dst.PS_value = static_cast<NUMERIC_TYPE*>(
			cuda::malloc_device(numPS*sizeof(NUMERIC_TYPE)));
	cuda::copy(dst.PS_value, src.PS_Val, numPS*sizeof(NUMERIC_TYPE));

	dst.PS_time_series = static_cast<TimeSeries**>(
			cuda::malloc_device(numPS*sizeof(TimeSeries*)));
	TimeSeries** time_series = new TimeSeries*[numPS];

	for (int ps_i=0; ps_i<numPS; ps_i++)
	{
		if (src.PS_Ident[ps_i] == HVAR3 || src.PS_Ident[ps_i] == QVAR5)
		{
			::TimeSeries* target = src.PS_TimeSeries[ps_i];

			for (std::vector<::TimeSeries>::size_type ts_i=0,
					e=src.allTimeSeries.size(); ts_i!=e; ts_i++)
			{
				::TimeSeries* candidate = &(src.allTimeSeries[ts_i]);
				if (candidate->value == target->value)
				{
					time_series[ps_i] = &(dst.all_time_series[ts_i]);
						//+ ts_i*sizeof(TimeSeries*);
				}
			}
		}
	}
	cuda::copy(dst.PS_time_series,
			time_series, numPS*sizeof(TimeSeries*));

	dst.PS_idx = static_cast<int*>(
			cuda::malloc_device(numPS*sizeof(int)));
	int* idx = new int[numPS];

	for (int ps_i=0; ps_i<numPS; ps_i++)
	{
		idx[ps_i] = src.ypi[ps_i] * pitch + src.xpi[ps_i] + offset;
	}

	cuda::copy(dst.PS_idx, idx, numPS*sizeof(int));
	delete[] time_series;
	delete[] idx;
}

void lis::cuda::Boundary::update_time_series
(
	NUMERIC_TYPE t
)
{
	cuda::update_time_series<<<1, 1>>>(t);
	update_time_varying_boundary_conditions<<<64, CUDA_BLOCK_SIZE>>>();
	update_time_varying_point_sources<<<1, CUDA_BLOCK_SIZE>>>();
}

void lis::cuda::Boundary::update_point_sources
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats
)
{
	cuda::update_point_sources<<<1, CUDA_BLOCK_SIZE>>>(H, DEM, mass_stats);
}

void lis::cuda::Boundary::update_H
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* Qx,
	NUMERIC_TYPE* Qy,
	MassStats* mass_stats,
	dim3 grid_size,
	int drycheck
)
{
	cuda::update_point_sources_Q << <1, CUDA_BLOCK_SIZE >> > (H, mass_stats);
	cuda::update_H << <grid_size, cuda::block_size >> > (H, Qx, Qy, drycheck);
	cuda::update_point_sources_H << <1, CUDA_BLOCK_SIZE >> > (H, DEM, mass_stats);
}								  
void lis::cuda::Boundary::drain_nodata_water
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats,
	dim3 grid_size
)
{
	cuda::drain_nodata_water<<<grid_size, CUDA_BLOCK_SIZE>>>(H, DEM,
			mass_stats);
}

void lis::cuda::Boundary::drain_nodata_waterACC
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	MassStats* mass_stats,
	dim3 grid_size
)
{
	cuda::drain_nodata_waterACC<<<grid_size, CUDA_BLOCK_SIZE>>>(H, DEM,
			mass_stats);
}

void lis::cuda::Boundary::free_device
(
	BoundaryConditions& boundaries
)
{
	for (int i=0; i<boundaries.time_series_count; i++)
	{
		cuda::free_device(boundaries.all_time_series[i].time);
		cuda::free_device(boundaries.all_time_series[i].value);
	}

	cuda::free_device(boundaries.BC_type);
	cuda::free_device(boundaries.BC_value);
	cuda::free_device(boundaries.BC_time_series);
	cuda::free_device(boundaries.PS_type);
	cuda::free_device(boundaries.PS_value);
	cuda::free_device(boundaries.PS_time_series);
	cuda::free_device(boundaries.PS_idx);
	cuda::free_device(boundaries.all_time_series);
}

