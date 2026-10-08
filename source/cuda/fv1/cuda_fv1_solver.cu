#include "cuda_fv1_solver.cuh"
#include "cuda_boundary.cuh"
#include "cuda_hll.cuh"
#include "cuda_solver.cuh"
#include <algorithm>

namespace lis
{
namespace cuda
{
namespace fv1
{

__global__ void zero_ghost_cells_north_south
(
	Flow U
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;

	for (int i=global_i; i<cuda::pitch; i+=blockDim.x*gridDim.x)
	{
		{
			int j=0;
			U.H[j*pitch + i] = C(0.0);
			U.HU[j*pitch + i] = C(0.0);
			U.HV[j*pitch + i] = C(0.0);
		}

		{
			int j=cuda::geometry.ysz+1;
			U.H[j*pitch + i] = C(0.0);
			U.HU[j*pitch + i] = C(0.0);
			U.HV[j*pitch + i] = C(0.0);
		}
	}
}

__global__ void zero_ghost_cells_east_west
(
	Flow U
)
{
	int global_j = blockIdx.x*blockDim.x + threadIdx.x;

	for (int j=global_j; j<cuda::geometry.ysz+2; j+=blockDim.x*gridDim.x)
	{
		{
			int i=0;
			U.H[j*pitch + i] = C(0.0);
			U.HU[j*pitch + i] = C(0.0);
			U.HV[j*pitch + i] = C(0.0);
		}

		{
			int i=cuda::geometry.xsz+1;
			U.H[j*pitch + i] = C(0.0);
			U.HU[j*pitch + i] = C(0.0);
			U.HV[j*pitch + i] = C(0.0);
		}
	}
}

__global__ void update_dt_per_element
(
	NUMERIC_TYPE* dt,
	Flow U
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j; j<cuda::geometry.ysz+2; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::pitch; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE H = U.H[j*cuda::pitch + i];

			if (H > cuda::solver_params.DepthThresh)
			{
				NUMERIC_TYPE HU = U.HU[j*cuda::pitch + i];
				NUMERIC_TYPE HV = U.HV[j*cuda::pitch + i];
				NUMERIC_TYPE U = HU/H;
				NUMERIC_TYPE V = HV/H;
				
				NUMERIC_TYPE dt_x = cuda::solver_params.cfl * cuda::geometry.dx
					/ (FABS(U)+SQRT(cuda::physical_params.g * H));
				NUMERIC_TYPE dt_y = cuda::solver_params.cfl * cuda::geometry.dy
					/ (FABS(V)+SQRT(cuda::physical_params.g * H));

				dt[j*cuda::pitch + i] = FMIN(dt_x, dt_y);
			}
			else
			{
				dt[j*cuda::pitch + i] = cuda::solver_params.max_dt;
			}
		}
	}
}

__global__
void update_uniform_rain_func
(
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE  rain_rate
)
{
	int global_i = blockIdx.x * blockDim.x + threadIdx.x;
	int global_j = blockIdx.y * blockDim.y + threadIdx.y;

	for (int j = global_j; j < cuda::geometry.ysz+2; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::pitch; i += blockDim.x * gridDim.x)
		{
			NUMERIC_TYPE cell_rain;

			cell_rain = rain_rate * cuda::dt;

			NUMERIC_TYPE Z = DEM[j * cuda::pitch + i];
			if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6) /* || Z < C(0.0) */) {

			}
			else {

				NUMERIC_TYPE& Hval = U.H[j * cuda::pitch + i];
				Hval += cell_rain;
			}
		}
	}
}

__global__ void update_ghost_cells
(
	Flow U,
	NUMERIC_TYPE* Zstar_x,
	NUMERIC_TYPE* Zstar_y
)
{
	for (int j=blockIdx.x*blockDim.x+threadIdx.x+1; j<cuda::geometry.ysz+1;
			j+=blockDim.x*gridDim.x)
	{
		// west
		{
			int i = 1;
			NUMERIC_TYPE Zstar = Zstar_x[j*cuda::pitch + i-1];
			FlowVector U_inside = U[j*cuda::pitch + i];

			FlowVector U_outside = Boundary::outside_x(U_inside, U_inside,
					Boundary::index_w(i, j), 1, Zstar);

			U.H[j*cuda::pitch + i-1] = U_outside.H;
			U.HU[j*cuda::pitch + i-1] = U_outside.HU;
			U.HV[j*cuda::pitch + i-1] = U_outside.HV;
		}

		// east
		{
			int i = cuda::geometry.xsz;
			NUMERIC_TYPE Zstar = Zstar_x[j*cuda::pitch + i];
			FlowVector U_inside = U[j*cuda::pitch + i];

			FlowVector U_outside = Boundary::outside_x(U_inside, U_inside,
					Boundary::index_e(i, j), -1, Zstar);

			U.H[j*cuda::pitch + i+1] = U_outside.H;
			U.HU[j*cuda::pitch + i+1] = U_outside.HU;
			U.HV[j*cuda::pitch + i+1] = U_outside.HV;
		}
	}

	for (int i=blockIdx.x*blockDim.x+threadIdx.x+1; i<cuda::geometry.xsz+1;
			i+=blockDim.x*gridDim.x)
	{
		// north
		{
			int j = 1;
			NUMERIC_TYPE Zstar = Zstar_y[(j-1)*cuda::pitch + i];
			FlowVector U_inside = U[j*cuda::pitch + i];

			FlowVector U_outside = Boundary::outside_y(U_inside, U_inside,
					Boundary::index_n(i, j), -1, Zstar);

			U.H[(j-1)*cuda::pitch + i] = U_outside.H;
			U.HU[(j-1)*cuda::pitch + i] = U_outside.HU;
			U.HV[(j-1)*cuda::pitch + i] = U_outside.HV;
		}

		// south
		{
			int j = cuda::geometry.ysz;
			NUMERIC_TYPE Zstar = Zstar_y[j*cuda::pitch + i];
			FlowVector U_inside = U[j*cuda::pitch + i];

			FlowVector U_outside = Boundary::outside_y(U_inside, U_inside,
					Boundary::index_s(i, j), 1, Zstar);

			U.H[(j+1)*cuda::pitch + i] = U_outside.H;
			U.HU[(j+1)*cuda::pitch + i] = U_outside.HU;
			U.HV[(j+1)*cuda::pitch + i] = U_outside.HV;
		}
	}
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
apply_friction
(
	Flow U,
	NUMERIC_TYPE* manning
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	for (int j=global_j+1; j<cuda::geometry.ysz+1; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i+1; i<cuda::geometry.xsz+1; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE H = U.H[j*cuda::pitch + i];
			NUMERIC_TYPE& HU = U.HU[j*cuda::pitch + i];
			NUMERIC_TYPE& HV = U.HV[j*cuda::pitch + i];

			if (H <= cuda::solver_params.DepthThresh) {
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

			NUMERIC_TYPE U = HU/H;
			NUMERIC_TYPE V = HV/H;
			if (FABS(U) <= cuda::solver_params.SpeedThresh
					&& FABS(V) <= cuda::solver_params.SpeedThresh)
			{
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

            NUMERIC_TYPE n = (manning == nullptr)
                ? cuda::physical_params.manning
				: manning[j*cuda::pitch + i];

			NUMERIC_TYPE Cf = cuda::physical_params.g * n * n /
				POW(H, C(1.0)/C(3.0));
			NUMERIC_TYPE speed = SQRT(U*U+V*V);

			NUMERIC_TYPE Sf_x = -Cf*U*speed;
			NUMERIC_TYPE Sf_y = -Cf*V*speed;
			NUMERIC_TYPE D_x = C(1.0) + cuda::dt*Cf / H
				* (C(2.0)*U*U+V*V)/speed;
			NUMERIC_TYPE D_y = C(1.0) + cuda::dt*Cf / H
				* (U*U+C(2.0)*V*V)/speed;

			HU += cuda::dt*Sf_x/D_x;
			HV += cuda::dt*Sf_y/D_y;
		}
	}
}

__device__ NUMERIC_TYPE bed_source_x
(
 	NUMERIC_TYPE Zstar_w,
	NUMERIC_TYPE Zstar_e,
	NUMERIC_TYPE Hstar_w,
	NUMERIC_TYPE Hstar_e,
	NUMERIC_TYPE ETA
)
{
	NUMERIC_TYPE Zdagger_w = Zstar_w - FMAX(C(0.0), -(ETA - Zstar_w));
	NUMERIC_TYPE Zdagger_e = Zstar_e - FMAX(C(0.0), -(ETA - Zstar_e));

	return -cuda::physical_params.g * C(0.5) * (Hstar_w + Hstar_e)
		* (Zdagger_e - Zdagger_w)/cuda::geometry.dx;
}

__device__ NUMERIC_TYPE bed_source_y
(
 	NUMERIC_TYPE Zstar_s,
	NUMERIC_TYPE Zstar_n,
	NUMERIC_TYPE Hstar_s,
	NUMERIC_TYPE Hstar_n,
	NUMERIC_TYPE ETA
)
{
	NUMERIC_TYPE Zdagger_s = Zstar_s - FMAX(C(0.0), -(ETA - Zstar_s));
	NUMERIC_TYPE Zdagger_n = Zstar_n - FMAX(C(0.0), -(ETA - Zstar_n));

	return -cuda::physical_params.g * C(0.5) * (Hstar_s + Hstar_n)
		* (Zdagger_n - Zdagger_s)/cuda::geometry.dy;
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
update_flow_variables_x
(
	Flow Uold,
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* Zstar_x,
	MassStats* mass_stats
)
{
	__shared__ FlowVector F[CUDA_BLOCK_SIZE_Y][CUDA_BLOCK_SIZE_X];
	__shared__ NUMERIC_TYPE Hstar[CUDA_BLOCK_SIZE_Y][CUDA_BLOCK_SIZE_X];

	int global_i = blockIdx.x*(blockDim.x-1) + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

	int max_i = ((cuda::geometry.xsz+blockDim.x-2)/(blockDim.x-1))*(blockDim.x-1);
	int max_j = ((cuda::geometry.ysz+blockDim.y-1)/blockDim.y)*blockDim.y;

	for (int j=global_j+1; j<=max_j; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<=max_i; i+=(blockDim.x-1)*gridDim.x)
		{
			// blocks overlap in the x direction
			// so don't start a new block at the right-hand edge
			if (threadIdx.x == 0 && i == max_i) continue;

			FlowVector U_neg;
			FlowVector Ustar_neg;
			FlowVector F_e;
			NUMERIC_TYPE Z_neg;
			NUMERIC_TYPE Zstar_e;

			// only update flow variables if (i, j) is inside the domain
			if (i <= cuda::geometry.xsz && j <= cuda::geometry.ysz)
			{
				Zstar_e = Zstar_x[j*cuda::pitch + i];

				Z_neg = DEM[j*cuda::pitch + i];
				U_neg = Uold[j*cuda::pitch + i];
				Ustar_neg = U_neg.star(Z_neg, Zstar_e);

				NUMERIC_TYPE Z_pos = DEM[j*cuda::pitch + i+1];
				FlowVector U_pos = Uold[j*cuda::pitch + i+1];
				FlowVector Ustar_pos = U_pos.star(Z_pos, Zstar_e);

				if (i == 0)
				{
					Ustar_pos = Boundary::inside_x(Ustar_neg, Ustar_pos,
							Ustar_pos, Boundary::index_w(i, j));
				}
				else if (i == cuda::geometry.xsz)
				{
					Ustar_neg = Boundary::inside_x(Ustar_pos, Ustar_neg,
							Ustar_neg, Boundary::index_e(i, j));
				}

				Hstar[threadIdx.y][threadIdx.x] = Ustar_pos.H;

				F_e = HLL::x(Ustar_neg, Ustar_pos);
				F[threadIdx.y][threadIdx.x] = F_e;
				update_mass_stats_x(mass_stats, F_e.H, i, j);

			}

			__syncthreads();

			if (threadIdx.x == 0) continue;

			if (i <= cuda::geometry.xsz && j <= cuda::geometry.ysz)
			{
				FlowVector& F_w = F[threadIdx.y][threadIdx.x-1];
				FlowVector U0 = Uold[j*cuda::pitch + i];

				NUMERIC_TYPE& H = U.H[j*cuda::pitch + i];
				NUMERIC_TYPE& HU = U.HU[j*cuda::pitch + i];
				NUMERIC_TYPE& HV = U.HV[j*cuda::pitch + i];

				NUMERIC_TYPE Zstar_w = Zstar_x[j*cuda::pitch + i-1];
				NUMERIC_TYPE Hstar_w = Hstar[threadIdx.y][threadIdx.x-1];
				NUMERIC_TYPE ETA = U_neg.H + Z_neg;

				H = U0.H - cuda::dt * (F_e.H - F_w.H)/cuda::geometry.dx;
				HU = U0.HU - cuda::dt * ((F_e.HU - F_w.HU)/cuda::geometry.dx
					- bed_source_x(Zstar_w, Zstar_e, Hstar_w, Ustar_neg.H, ETA));
				HV = U0.HV - cuda::dt * (F_e.HV - F_w.HV)/cuda::geometry.dx;

			}
		}
	}
}

__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
update_flow_variables_y
(
	Flow Uold,
	Flow Uint,
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* Zstar_y,
	MassStats* mass_stats
)
{
        __shared__ FlowVector F[CUDA_BLOCK_SIZE_Y][CUDA_BLOCK_SIZE_X];
        __shared__ NUMERIC_TYPE Hstar[CUDA_BLOCK_SIZE_Y][CUDA_BLOCK_SIZE_X];

        int global_i = blockIdx.x*blockDim.x + threadIdx.x;
        int global_j = blockIdx.y*(blockDim.y-1) + threadIdx.y;

	int max_i = ((cuda::geometry.xsz+blockDim.x-1)/blockDim.x)*blockDim.x;
	int max_j = ((cuda::geometry.ysz+blockDim.y-2)/(blockDim.y-1))*(blockDim.y-1);

	for (int j=global_j; j<=max_j; j+=(blockDim.y-1)*gridDim.y)
	{
		if (threadIdx.y == 0 && j == max_j) continue;

		for (int i=global_i+1; i<=max_i; i+=blockDim.x*gridDim.x)
		{
			FlowVector U_pos;
			FlowVector Ustar_pos;
			FlowVector F_s;
			NUMERIC_TYPE Z_pos;
			NUMERIC_TYPE Zstar_s;

			if (i <= cuda::geometry.xsz && j <= cuda::geometry.ysz)
			{
				Zstar_s = Zstar_y[j*cuda::pitch + i];

				NUMERIC_TYPE Z_neg = DEM[(j+1)*cuda::pitch + i];
				FlowVector U_neg = Uold[(j+1)*cuda::pitch + i];
				FlowVector Ustar_neg = U_neg.star(Z_neg, Zstar_s);

				Z_pos = DEM[j*cuda::pitch + i];
				U_pos = Uold[j*cuda::pitch + i];
				Ustar_pos = U_pos.star(Z_pos, Zstar_s);

				// this hack appears necessary to prevent some unknown
				// compiler optimisation from causing a bug seen in
				// EA test 4
				//Ustar_pos.H = Ustar_pos.H*C(2.0)/C(2.0);

				if (j == 0)
				{       
					Ustar_neg = Boundary::inside_y(Ustar_pos, Ustar_neg,
							Ustar_neg, Boundary::index_n(i, j));
				}
				else if (j == cuda::geometry.ysz)
				{
					Ustar_pos = Boundary::inside_y(Ustar_neg, Ustar_pos,
							Ustar_pos, Boundary::index_s(i, j));
				}

				Hstar[threadIdx.y][threadIdx.x] = Ustar_neg.H;

				F_s = HLL::y(Ustar_neg, Ustar_pos);
				F[threadIdx.y][threadIdx.x] = F_s;
				update_mass_stats_y(mass_stats, F_s.H, i, j);
			}

                        __syncthreads();

                        if (threadIdx.y == 0) continue;

			if (i <= cuda::geometry.xsz && j <= cuda::geometry.ysz)
			{
				FlowVector& F_n = F[threadIdx.y-1][threadIdx.x];
				FlowVector U0 = Uint[j*cuda::pitch + i];

				NUMERIC_TYPE& H = U.H[j*cuda::pitch + i];
				NUMERIC_TYPE& HU = U.HU[j*cuda::pitch + i];
				NUMERIC_TYPE& HV = U.HV[j*cuda::pitch + i];

				NUMERIC_TYPE Zstar_n = Zstar_y[(j-1)*cuda::pitch + i];
				NUMERIC_TYPE Hstar_n = Hstar[threadIdx.y-1][threadIdx.x];
				NUMERIC_TYPE ETA = U_pos.H + Z_pos;

				H = U0.H - cuda::dt * (F_n.H - F_s.H)/cuda::geometry.dy;
				HU = U0.HU - cuda::dt * (F_n.HU - F_s.HU)/cuda::geometry.dy;
				HV = U0.HV - cuda::dt * ((F_n.HV - F_s.HV)/cuda::geometry.dy
					- bed_source_y(Zstar_s, Zstar_n, Ustar_pos.H, Hstar_n, ETA));
			
			}
        }
    }
}

}
}
}

lis::cuda::fv1::Solver::Solver
(
	Flow& U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* Zstar_x,
	NUMERIC_TYPE* Zstar_y,
	NUMERIC_TYPE* manning,
	Geometry& geometry,
	PhysicalParams& physical_params,
	dim3 grid_size
)
:
Uold(U1),
U(U2),
DEM(DEM),
Zstar_x(Zstar_x),
Zstar_y(Zstar_y),
manning(manning),
grid_size(grid_size)
{
	Flow::allocate_device(U1, geometry);
	Flow::allocate_device(U2, geometry);
	Flow::allocate_device(Ux, geometry);

	Flow::copy(U1, U, geometry);

	friction = (physical_params.manning > C(0.0)) || (manning != nullptr);
}

void lis::cuda::fv1::Solver::zero_ghost_cells()
{
	zero_ghost_cells_north_south<<<1, CUDA_BLOCK_SIZE>>>(U);
	zero_ghost_cells_east_west<<<1, CUDA_BLOCK_SIZE>>>(U);
}

void lis::cuda::fv1::Solver::update_ghost_cells()
{
	lis::cuda::fv1::update_ghost_cells<<<64, CUDA_BLOCK_SIZE>>>(
			Uold, Zstar_x, Zstar_y);
}

lis::cuda::fv1::Flow& lis::cuda::fv1::Solver::update_flow_variables
(
	MassStats* mass_stats
)
{
	if (friction)
	{
		apply_friction<<<grid_size, cuda::block_size>>>(Uold, manning);
		update_ghost_cells();
	}

	update_flow_variables_x<<<grid_size, cuda::block_size>>>
			(Uold, Ux, DEM, Zstar_x, mass_stats);

	update_flow_variables_y<<<grid_size, cuda::block_size>>>
			(Uold, Ux, U, DEM, Zstar_y, mass_stats);

	std::swap(Uold, U);

	return Uold;
}

lis::cuda::fv1::Flow& lis::cuda::fv1::Solver::d_U()
{
	return Uold;
}

void lis::cuda::fv1::Solver::update_dt_per_element
(
	NUMERIC_TYPE* dt_field
) const
{
	lis::cuda::fv1::update_dt_per_element<<<grid_size, cuda::block_size>>>(
			dt_field, Uold);
}



void lis::cuda::fv1::Solver::FloodplainQ()
{

  return;	
}

void lis::cuda::fv1::Solver::zero_thin_depth_slopes()
{

  return;	
}

void lis::cuda::fv1::Solver::updateMaxFieldACC
(
	NUMERIC_TYPE t
)
{
	return;
}

void lis::cuda::fv1::Solver::update_uniform_rain
(
	NUMERIC_TYPE rain_rate
)
{
	update_uniform_rain_func << <grid_size, cuda::block_size >> > (
		Uold, DEM, rain_rate);
}


lis::cuda::fv1::Solver::~Solver()
{
	Flow::free_device(U1);
	Flow::free_device(U2);
	Flow::free_device(Ux);
}
