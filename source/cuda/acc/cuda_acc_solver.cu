#include "cuda_acc_solver.cuh"
#include "cuda_boundary.cuh"
#include "cuda_hll.cuh"
#include "cuda_solver.cuh"
#include <algorithm>

__managed__ NUMERIC_TYPE lis::cuda::acc::maxH;

namespace lis
{
namespace cuda
{
namespace acc
{

//	// Calculate timestep for acceleration version based on 2D shallow water CFL condition
//__device__ void CalcT 
//(
//	Flow U
//	
//)
//{
//	NUMERIC_TYPE cfl=C(0.7), g;
//	NUMERIC_TYPE MH=C(0.0);
//	NUMERIC_TYPE locT;
//
//	g=cuda::physical_params.g;
//	cfl=cuda::solver_params.cfl;
//	
//	
//	MH=cuda::acc::maxH; 
//	
//	// Apply local 2D CFL condition and use global minimum timestep to ensure global stability
//	if(MH >	cuda::solver_params.DepthThresh) //Don't apply where h=0 as h appears in equation denominator
//	{
//		// Time step control for stability from actual equations (MSH, implemented by TJF)
//		locT=cfl*cuda::geometry.dx/(SQRT(g*MH));
//		cuda::dt=getmin(cuda::dt,locT);   
//	}
//	else // Set to initial timestep if h=0
//	{
//		cuda::dt=cuda::solver_params.InitTstep;
//	}
//	
//return;
//}



__global__ void zero_ghost_cells_north_south
(
	Flow U
)
{
	int global_i = blockIdx.x * blockDim.x + threadIdx.x;

	for (int i = global_i; i < cuda::pitch; i += blockDim.x * gridDim.x)
	{
		{
			int j = 0;
			U.H[j * pitch + i] = C(0.0);

		}

		{
			int j = cuda::geometry.ysz + 1;
			U.H[j * pitch + i] = C(0.0);

		}
	}
}

__global__ void zero_ghost_cells_east_west
(
	Flow U
)
{
	int global_j = blockIdx.x * blockDim.x + threadIdx.x;

	for (int j = global_j; j < cuda::geometry.ysz + 2; j += blockDim.x * gridDim.x)
	{
		{
			int i = 0;
			U.H[j * pitch + i] = C(0.0);

		}

		{
			int i = cuda::geometry.xsz + 1;
			U.H[j * pitch + i] = C(0.0);

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

	for (int j=global_j; j<cuda::geometry.ysz; j+=blockDim.y*gridDim.y)
	{
		for (int i=global_i; i<cuda::pitch; i+=blockDim.x*gridDim.x)
		{
			NUMERIC_TYPE H = U.H[j*cuda::pitch + i];

			if (H > cuda::solver_params.DepthThresh)
			{


				NUMERIC_TYPE dt_ = cuda::solver_params.cfl * cuda::geometry.dx
					/ SQRT(cuda::physical_params.g * H);

				if (FABS(dt[j * cuda::pitch + i] - C(0.0)) < C(1e-6)) {
					dt[j * cuda::pitch + i] = cuda::solver_params.max_dt;
				}
				else {
					dt[j * cuda::pitch + i] = getmin(dt_, dt[j * cuda::pitch + i]);
				}


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

	for (int j = global_j; j < cuda::geometry.ysz; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::geometry.xsz; i += blockDim.x * gridDim.x)
		{
			NUMERIC_TYPE cell_rain;

			cell_rain = rain_rate * cuda::dt;

			NUMERIC_TYPE Z = DEM[j * cuda::geometry.xsz + i];
			if (FABS(Z - cuda::solver_params.nodata_elevation) < C(1e-6) /* || Z < C(0.0) */) {

			}
			else {

				NUMERIC_TYPE& Hval = U.H[j * cuda::geometry.xsz + i];
				Hval += cell_rain;
			}
		}
	}
}

__global__ void update_ghost_cells
(
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* manning
)
{
	int p0, p1, sign, edge, dir;
	NUMERIC_TYPE h0, h1, z0, z1, hflow, fn, dh, Sf, g, q0;
	NUMERIC_TYPE qoldptr;
	
	g = cuda::physical_params.g;
    dt = cuda::dt;
	
	for (int j=blockIdx.x*blockDim.x+threadIdx.x; j<cuda::geometry.ysz;
			j+=blockDim.x*gridDim.x)
	{
		{
		    U.Qx[j*(cuda::geometry.xsz + 1)] = C(0.0);  //zero West boundary
		    U.Qx[cuda::geometry.xsz + j*(cuda::geometry.xsz + 1)] = C(0.0); //zero East boundary			
		}	
	}
	
	for (int i=blockIdx.x*blockDim.x+threadIdx.x; i<cuda::geometry.xsz;
			i+=blockDim.x*gridDim.x)
	{
		{
		    U.Qy[i] = C(0.0);  //zero North boundary
		    U.Qy[i + cuda::geometry.ysz*(cuda::geometry.xsz + 1)] = C(0.0); //zero South boundary
		}
	}	
					
	for (int j=blockIdx.x*blockDim.x+threadIdx.x; j<cuda::geometry.ysz;
			j+=blockDim.x*gridDim.x)
	{
		// west
		{
			int i = 1;
			int bc_i = Boundary::index_w_ACC(i, j);
            p0 = j * cuda::geometry.xsz;
			p1 = p0 + 1;

            sign = -1;
			NUMERIC_TYPE& qptr = U.Qx[j * (cuda::geometry.xsz + 1)];
            qoldptr = U.Qxold[j * (cuda::geometry.xsz + 1)];
            q0 = qoldptr; // Set old value of Q for acceleration version
            edge = 4;
			dir = 2;
			
			if (cuda::boundaries.BC_type[bc_i] == FREE1 && U.H[p0] > cuda::solver_params.DepthThresh)
			{	// FREE boundary
				hflow = U.H[p0];
				h0 = U.H[p0];
				h1 = U.H[p1];
				z0 = DEM[p0];
				z1 = DEM[p1];
				if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;				
                
                if (cuda::boundaries.BC_value[bc_i] < C(-0.999)) 
                {
                	dh = z0 + h0 - z1 - h1;
                	 
                	Sf = -dh / cuda::geometry.dx;
                }
                else // use floodplain slope
                {
                	Sf = cuda::boundaries.BC_value[bc_i];
                }
                qoldptr = sign*(FABS(q0) + FABS(g*cuda::dt*hflow*Sf)) / (C(1.0) + g*cuda::dt*hflow*fn*fn*FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
                qptr = qoldptr*cuda::geometry.dx; // potntially unnessasary?? check for repetition in UpdateQs (jcn)

				if (dir == 1 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 1 && sign == 1 && qptr < 0) qptr = C(0.0);
				if (dir == 2 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 2 && sign == 1 && qptr < 0) qptr = C(0.0);
			}				
			
			
            switch (cuda::boundaries.BC_type[bc_i])
            {
            case HFIX2:
            case HVAR3:
                if (U.H[p0] > cuda::solver_params.DepthThresh)
                {
                	h0 = cuda::boundaries.BC_value[bc_i];
                	hflow = getmax(U.H[p0], h0 - DEM[p0]);
                	h1 = U.H[p0];
                	z1 = DEM[p0];
                
                	if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;
                
                	//if (!(U.ChanMask[p0] != -1)){
                		qptr = C(0.0);
                		dh = h0 - z1 - h1;
                		// change slops direction depending on the edge
                		if (edge == 1 || edge == 4) Sf = -dh / cuda::geometry.dx;
                		else Sf = dh / cuda::geometry.dx;
                		qoldptr = (q0 - g*cuda::dt*hflow*Sf) / (C(1.0) + g*cuda::dt*hflow*fn*fn*FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
                		qptr = qoldptr*cuda::geometry.dx;
                	//}
                }
                break; 
            // QVAR boundary
            case QVAR5:
            	{
            		h0 = cuda::boundaries.BC_value[bc_i];
            		qptr = -sign*h0*cuda::geometry.dx;
            		
            	}
            	break;
            default:
            	break;
            }
		}

		// east
		{
			int i = cuda::geometry.xsz;
			int bc_i = Boundary::index_e_ACC(i, j);
            p0 = cuda::geometry.xsz - 1 + j * cuda::geometry.xsz;
			p1 = p0 - 1;
            sign = 1;
			NUMERIC_TYPE& qptr = U.Qx[cuda::geometry.xsz + j * (cuda::geometry.xsz + 1)];
            qoldptr = U.Qxold[cuda::geometry.xsz + j * (cuda::geometry.xsz + 1)];
            
            q0 = qoldptr; // Set old value of Q for acceleration version
            edge = 2;
			dir = 2;
			
			if (cuda::boundaries.BC_type[bc_i] == FREE1 && U.H[p0] > cuda::solver_params.DepthThresh)
			{	// FREE boundary
				hflow = U.H[p0];
				h0 = U.H[p0];
				h1 = U.H[p1];
				z0 = DEM[p0];
				z1 = DEM[p1];
				if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;

				if (cuda::boundaries.BC_value[bc_i] < C(-0.999))
				{
					dh = z0 + h0 - z1 - h1;

					Sf = -dh / cuda::geometry.dx;
				}
				else // use floodplain slope
				{
					Sf = cuda::boundaries.BC_value[bc_i];
				}
				qoldptr = sign * (FABS(q0) + FABS(g * cuda::dt * hflow * Sf)) / (C(1.0) + g * cuda::dt * hflow * fn * fn * FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
				qptr = qoldptr * cuda::geometry.dx; // potntially unnessasary?? check for repetition in UpdateQs (jcn)	

				if (dir == 1 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 1 && sign == 1 && qptr < 0) qptr = C(0.0);
				if (dir == 2 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 2 && sign == 1 && qptr < 0) qptr = C(0.0);
			}
			
            switch (cuda::boundaries.BC_type[bc_i])
            {
            case HFIX2:
            case HVAR3:
                if (U.H[p0] > cuda::solver_params.DepthThresh)
                {
                	h0 = cuda::boundaries.BC_value[bc_i];
                	hflow = getmax(U.H[p0], h0 - DEM[p0]);
                	h1 = U.H[p0];
                	z1 = DEM[p0];
                
                	if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;
                
                	//if (!(U.ChanMask[p0] != -1)){
                		qptr = C(0.0);
                		dh = h0 - z1 - h1;
                		// change slops direction depending on the edge
                		if (edge == 1 || edge == 4) Sf = -dh / cuda::geometry.dx;
                		else Sf = dh / cuda::geometry.dx;
                		qoldptr = (q0 - g*cuda::dt*hflow*Sf) / (C(1.) + g*cuda::dt*hflow*fn*fn*FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
                		qptr = qoldptr*cuda::geometry.dx;
                	//}
                }
                break; 
            // QVAR boundary
            case QVAR5:
            	{
            		h0 = cuda::boundaries.BC_value[bc_i];
            		qptr = -sign*h0*cuda::geometry.dx;
            		
            	}
            	break;
            default:
            	break;
            }			
		}	
	}			
			
	for (int i=blockIdx.x*blockDim.x+threadIdx.x; i<cuda::geometry.xsz;
			i+=blockDim.x*gridDim.x)
	{
		// north
		{
			int j = 1;
			int bc_i = Boundary::index_n_ACC(i, j);
			p0 = i; 
			p1 = i + cuda::geometry.xsz;
		
			sign = -1;
			NUMERIC_TYPE& qptr = U.Qy[i]; 
			qoldptr = U.Qyold[i]; 
			q0 = qoldptr; 
			edge = 1;
			dir = 1;

			if (cuda::boundaries.BC_type[bc_i] == FREE1 && U.H[p0] > cuda::solver_params.DepthThresh)
			{	// FREE boundary
				hflow = U.H[p0];
				h0 = U.H[p0];
				h1 = U.H[p1];
				z0 = DEM[p0];
				z1 = DEM[p1];
				if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;

				if (cuda::boundaries.BC_value[bc_i] < C(-0.999))
				{
					dh = z0 + h0 - z1 - h1;

					Sf = -dh / cuda::geometry.dx;
				}
				else // use floodplain slope
				{
					Sf = cuda::boundaries.BC_value[bc_i];
				}
				qoldptr = sign * (FABS(q0) + FABS(g * cuda::dt * hflow * Sf)) / (C(1.0) + g * cuda::dt * hflow * fn * fn * FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
				qptr = qoldptr * cuda::geometry.dx; // potntially unnessasary?? check for repetition in UpdateQs (jcn)	

				if (dir == 1 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 1 && sign == 1 && qptr < 0) qptr = C(0.0);
				if (dir == 2 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 2 && sign == 1 && qptr < 0) qptr = C(0.0);
			}
			
            switch (cuda::boundaries.BC_type[bc_i])
            {
            case HFIX2:
            case HVAR3:
                if (U.H[p0] > cuda::solver_params.DepthThresh)
                {
                	h0 = cuda::boundaries.BC_value[bc_i];
                	hflow = getmax(U.H[p0], h0 - DEM[p0]);
                	h1 = U.H[p0];
                	z1 = DEM[p0];
                
                	if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;
                
                	//if (!(U.ChanMask[p0] != -1)){
                		qptr = C(0.0);
                		dh = h0 - z1 - h1;
                		// change slops direction depending on the edge
                		if (edge == 1 || edge == 4) Sf = -dh / cuda::geometry.dx;
                		else Sf = dh / cuda::geometry.dx;
                		qoldptr = (q0 - g*cuda::dt*hflow*Sf) / (C(1.0) + g*cuda::dt*hflow*fn*fn*FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
                		qptr = qoldptr*cuda::geometry.dx;
                	//}
                }
                break; 
            // QVAR boundary
            case QVAR5:
            	{
            		h0 = cuda::boundaries.BC_value[bc_i];
            		qptr = -sign*h0*cuda::geometry.dx;
            		
            	}
            	break;
            default:
            	break;
            }			
		}

		// south
		{
			int j = cuda::geometry.ysz;
			int bc_i = Boundary::index_s_ACC(i, j);
            p0 = cuda::geometry.xsz * (cuda::geometry.ysz - 1) + i;
			p1 = cuda::geometry.xsz * (cuda::geometry.ysz - 2) + i;
        
            sign = 1;
			NUMERIC_TYPE& qptr = U.Qy[(cuda::geometry.xsz + 1) * cuda::geometry.ysz + i];
            qoldptr = U.Qyold[(cuda::geometry.xsz + 1) * cuda::geometry.ysz + i];
            q0 = qoldptr; // Set old value of Q for acceleration version
            edge = 3;
			dir = 1;
			
			if (cuda::boundaries.BC_type[bc_i] == FREE1 && U.H[p0] > cuda::solver_params.DepthThresh)
			{	// FREE boundary
				hflow = U.H[p0];
				h0 = U.H[p0];
				h1 = U.H[p1];
				z0 = DEM[p0];
				z1 = DEM[p1];
				if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;

				if (cuda::boundaries.BC_value[bc_i] < C(-0.999))
				{
					dh = z0 + h0 - z1 - h1;

					Sf = -dh / cuda::geometry.dx;
				}
				else // use floodplain slope
				{
					Sf = cuda::boundaries.BC_value[bc_i];
				}
				qoldptr = sign * (FABS(q0) + FABS(g * cuda::dt * hflow * Sf)) / (C(1.0) + g * cuda::dt * hflow * fn * fn * FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
				qptr = qoldptr * cuda::geometry.dx; // potntially unnessasary?? check for repetition in UpdateQs (jcn)	

				//if (qptr < 0.5) {
				//	printf("qptr %lf\n", qptr);
				//	printf("q0 %lf\n", q0);
				//	printf("qoldptr %lf\n", qoldptr);
				//	printf("dx %lf\n", cuda::geometry.dx);
				//	printf("hflow %lf\n", hflow);
				//	printf("fn %lf\n", fn);
				//	printf("g %lf\n", g);
				//	printf("h0 %lf\n", h0);
				//	printf("h1 %lf\n", h1);
				//	printf("z0 %lf\n", z0);
				//	printf("z1 %lf\n", z1);
				//	printf("Sf %lf\n", Sf);
				//}

				if (dir == 1 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 1 && sign == 1 && qptr < 0) qptr = C(0.0);
				if (dir == 2 && sign == -1 && qptr > 0) qptr = C(0.0);
				if (dir == 2 && sign == 1 && qptr < 0) qptr = C(0.0);
			}
			
            switch (cuda::boundaries.BC_type[bc_i])
            {
            case HFIX2:
            case HVAR3:
                if (U.H[p0] > cuda::solver_params.DepthThresh)
                {
                	h0 = cuda::boundaries.BC_value[bc_i];
                	hflow = getmax(U.H[p0], h0 - DEM[p0]);
                	h1 = U.H[p0];
                	z1 = DEM[p0];
                
                	if (manning != NULL) fn = manning[p0]; else fn = cuda::physical_params.manning;
                
                	//if (!(U.ChanMask[p0] != -1)){
                		qptr = C(0.0);
                		dh = h0 - z1 - h1;
                		// change slops direction depending on the edge
                		if (edge == 1 || edge == 4) Sf = -dh / cuda::geometry.dx;
                		else Sf = dh / cuda::geometry.dx;
                		qoldptr = (q0 - g*cuda::dt*hflow*Sf) / (C(1.) + g*cuda::dt*hflow*fn*fn*FABS(q0) / (POW(hflow, (C(10.) / C(3.)))));
                		qptr = qoldptr*cuda::geometry.dx;
                	//}
                }
                break; 
            // QVAR boundary
            case QVAR5:
            	{
            		h0 = cuda::boundaries.BC_value[bc_i];
            		qptr = -sign*h0*cuda::geometry.dx;
            		
            	}
            	break;
            default:
            	break;
            }			
		}
	}
}



__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
update_max_field_ACC
(
	Flow Uold,
	NUMERIC_TYPE time
)
{
	int global_i = blockIdx.x * blockDim.x + threadIdx.x;
	int global_j = blockIdx.y * blockDim.y + threadIdx.y;

	for (int j = global_j; j < cuda::geometry.ysz; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::geometry.xsz; i += blockDim.x * gridDim.x)
		{
//			NUMERIC_TYPE& maxHvalue = maxH[j * cuda::pitch + i];
//			maxHvalue = FMAX(maxHvalue, H[j * cuda::pitch + i]);

			int pqptr = i + j * cuda::pitch;

			if (Uold.H[pqptr] > cuda::solver_params.DepthThresh)
			{
				if (Uold.initHtm[pqptr] == NULLVAL)
					Uold.initHtm[pqptr] = time / C(3600.0);
				
				Uold.totalHtm[pqptr] += cuda::dt / C(3600.0);
				// Update maximum water depths, and time of maximum (in hours) TDF updated only track maxH when h > DepthThresh
				if (Uold.H[pqptr] > Uold.maxH[pqptr])
				{
					Uold.maxH[pqptr] = Uold.H[pqptr];
					Uold.maxHtm[pqptr] = time / C(3600.0);
				}
			}
		

		}
	}
}





__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
update_flow_variables_x
(
	Flow Uold
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

    for (int j = global_j; j < cuda::geometry.ysz; j += blockDim.y*gridDim.y)
    {
        for (int i = global_i; i < cuda::geometry.xsz + 1; i += blockDim.x*gridDim.x)
        {		
			
		  int pqptr=i+j*(cuda::geometry.xsz+1);
		  // Divide Q by dx to return m2/s
		  Uold.Qxold[pqptr]=Uold.Qx[pqptr]*(C(1.0)/cuda::geometry.dx);			
			
		}
	}
}






__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
update_flow_variables_y
(
	Flow Uold
)
{
	int global_i = blockIdx.x*blockDim.x + threadIdx.x;
	int global_j = blockIdx.y*blockDim.y + threadIdx.y;

    for (int j = global_j; j < cuda::geometry.ysz + 1; j += blockDim.y*gridDim.y)
    {
        for (int i = global_i; i < cuda::geometry.xsz; i += blockDim.x*gridDim.x)
        {		
			
		  int pqptr=i+j*(cuda::geometry.xsz+1);
		  // Divide Q by dx to return m2/s
		  Uold.Qyold[pqptr]=Uold.Qy[pqptr]*(C(1.0)/cuda::geometry.dx);			
			
		}
	}
}


__device__ int MaskTest(int m1, int m2)
{
	if (m1 == -1 && m2 == -1) return(1);
	if (m1 == -1 && m2 > 0) return(1);
	if (m1 > 0 && m2 == -1) return(1);
	return(0);
}



__device__ NUMERIC_TYPE CalcFPQxAcc_
(
	int i,
	int j,
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* manning,
	bool friction
)
{

	NUMERIC_TYPE z0, z1, h0, h1, Sf = C(0.0), hflow, fn, dh = C(0.0), Q, q0, g, dt;
	int p0, p1, pq0;
	int pqy1, pqy2, pqy3, pqy4;
	NUMERIC_TYPE qy_avg, qvect;
	NUMERIC_TYPE qup, qdown, y0, y1;

	p0 = j * cuda::geometry.xsz + i;
	p1 = j * cuda::geometry.xsz + i + 1;

	pq0 = j * (cuda::geometry.xsz + 1) + i + 1;

	z0 = DEM[p0];
	z1 = DEM[p1];
	h0 = U.H[p0];
	h1 = U.H[p1];
	q0 = U.Qxold[pq0]; // in m2/s
	g = cuda::physical_params.g;
	dt = cuda::dt;

	if (friction) {
		pqy1 = j * (cuda::geometry.xsz + 1) + i;
		pqy2 = j * (cuda::geometry.xsz + 1) + (i + 1);
		pqy3 = (j + 1) * (cuda::geometry.xsz + 1) + i;
		pqy4 = (j + 1) * (cuda::geometry.xsz + 1) + (i + 1);
		qy_avg = (U.Qyold[pqy1] + U.Qyold[pqy2] + U.Qyold[pqy3] + U.Qyold[pqy4]) / C(4.0);
		qvect = SQRT(q0 * q0 + qy_avg * qy_avg);
	}
	else {
		qvect = q0;
	}


	if (manning != NULL) fn = C(0.5) * (manning[p0] + manning[p1]);
	else fn = cuda::physical_params.manning;
	y0 = z0 + h0;
	y1 = z1 + h1;

	dh = y0 - y1;

	Sf = -dh / cuda::geometry.dx;
	hflow = getmax(z0 + h0, z1 + h1) - getmax(z0, z1);
	hflow = getmax(hflow, C(0.0));
	hflow = getmin(hflow, cuda::solver_params.MaxHflow);

	qup = U.Qxold[pq0 - 1];
	qdown = U.Qxold[pq0 + 1];

	if (i == 0 && qup == 0)
	{
		qup = q0;
	}
	if (i == cuda::geometry.xsz - 2 && qdown == 0)
	{
		qdown = q0;
	}

	if ( /*MaskTest(U.ChanMask[p0], U.ChanMask[p1]) && */ hflow > cuda::solver_params.DepthThresh)
	{
//		if (cuda::solver_params.routing == OFF || FABS(Sf) < cuda::solver_params.RouteSfThresh) // CCS If routing scheme is enabled, only calc Q where Sf is below threshold value
//		{
			//Q-centred scheme
			Q = ((cuda::solver_params.theta * q0 + C(0.5) * (C(1.0) - cuda::solver_params.theta) * (qup + qdown)) - (g * dt * hflow * Sf)) / (C(1.0) + g * dt * hflow * fn * fn * FABS(qvect) / (POW(hflow, (C(10.0) / C(3.0))))) * cuda::geometry.dx;

			//Correction for situations where centred scheme introduces mass errors
			//if(copysign(1,Q)!= copysign(1,dh)){
			if (Q * dh < C(0.0)) { // version of line above that will compile on windows machine
				//Semi-implicit scheme of Bates et al (2010)
				Q = (q0 - (g * dt * hflow * Sf)) / (C(1.0) + g * dt * hflow * fn * fn * FABS(qvect) / (POW(hflow, (C(10.0) / C(3.0))))) * cuda::geometry.dx;
			}
//		}
//		else Q = C(0.0);
	}
	else Q = C(0.0);
	// option to save V's
	if (cuda::solver_params.voutput == ON)
	{
		if (Q != C(0.0))
		{
			U.Vx[pq0] = Q / cuda::geometry.dx / hflow;
		}
		else U.Vx[pq0] = C(0.0);
	}
	return(Q);
}





__device__ NUMERIC_TYPE CalcFPQyAcc_
(
	int i,
	int j,
	Flow U,
	NUMERIC_TYPE* DEM,
	NUMERIC_TYPE* manning,
	bool friction
)
{

	NUMERIC_TYPE z0, z1, h0, h1, Sf = C(0.0), hflow, fn, dh = C(0.0), Q, q0, g, dt;
	int p0, p1, pq0;
	int pqx1, pqx2, pqx3, pqx4;
	NUMERIC_TYPE qx_avg, qvect;
	NUMERIC_TYPE qup, qdown, y0, y1;

	p0 = j * cuda::geometry.xsz + i;
	p1 = (j + 1) * cuda::geometry.xsz + i;

	pq0 = (j + 1) * (cuda::geometry.xsz + 1) + i;

	z0 = DEM[p0];
	z1 = DEM[p1];
	h0 = U.H[p0];
	h1 = U.H[p1];
	q0 = U.Qyold[pq0]; // in m2/s
	g = cuda::physical_params.g;
	dt = cuda::dt;

	if (friction) {
		pqx1 = j * (cuda::geometry.xsz + 1) + i;
		pqx2 = j * (cuda::geometry.xsz + 1) + (i + 1);
		pqx3 = (j + 1) * (cuda::geometry.xsz + 1) + i;
		pqx4 = (j + 1) * (cuda::geometry.xsz + 1) + (i + 1);
		qx_avg = (U.Qxold[pqx1] + U.Qxold[pqx2] + U.Qxold[pqx3] + U.Qxold[pqx4]) / C(4.0);
		qvect = SQRT(q0 * q0 + qx_avg * qx_avg);
	}
	else {
		qvect = q0;
	}


	if (manning != NULL) fn = C(0.5) * (manning[p0] + manning[p1]);
	else fn = cuda::physical_params.manning;
	y0 = z0 + h0;
	y1 = z1 + h1;

	dh = y0 - y1;

	Sf = -dh / cuda::geometry.dx;
	hflow = getmax(z0 + h0, z1 + h1) - getmax(z0, z1);
	hflow = getmax(hflow, C(0.0));
	hflow = getmin(hflow, cuda::solver_params.MaxHflow);

	qup = U.Qyold[j * (cuda::geometry.xsz + 1) + i];
	qdown = U.Qyold[(j + 2) * (cuda::geometry.xsz + 1) + i];

	if (j == 0 && qup == 0)
	{
		qup = q0;
	}
	if (j == cuda::geometry.ysz - 2 && qdown == 0)
	{
		qdown = q0;
	}

	if (/*MaskTest(U.ChanMask[p0], U.ChanMask[p1]) &&*/ hflow > cuda::solver_params.DepthThresh)
	{
//		if (cuda::solver_params.routing == OFF || FABS(Sf) < cuda::solver_params.RouteSfThresh) // CCS If routing scheme is enabled, only calc Q where Sf is below threshold value
//		{
			//Q-centred scheme
			Q = ((cuda::solver_params.theta * q0 + C(0.5) * (C(1.0) - cuda::solver_params.theta) * (qup + qdown)) - (g * dt * hflow * Sf)) / (C(1.0) + g * dt * hflow * fn * fn * FABS(qvect) / (POW(hflow, (C(10.0) / C(3.0))))) * cuda::geometry.dx;

			//Correction for situations where centred scheme introduces mass errors
			//if(copysign(1,Q)!= copysign(1,dh)){
			if (Q * dh < C(0.0)) { // version of line above that will compile on windows machine
				//Semi-implicit scheme of Bates et al (2010)
				Q = (q0 - (g * dt * hflow * Sf)) / (C(1.0) + g * dt * hflow * fn * fn * FABS(qvect) / (POW(hflow, (C(10.0) / C(3.0))))) * cuda::geometry.dx;
			}
//		}
//		else Q = C(0.0);
	}
	else Q = C(0.0);
	// option to save V's
	if (cuda::solver_params.voutput == ON)
	{
		if (Q != C(0.0))
		{
			U.Vy[pq0] = Q / cuda::geometry.dx / hflow;
		}
		else U.Vy[pq0] = C(0.0);
	}
	return(Q);
}




__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
FloodplainQx
(
	
	Flow U,
	NUMERIC_TYPE* DEM,
	
	NUMERIC_TYPE* manning,
	
	bool friction
)
{

	int global_i = blockIdx.x * blockDim.x + threadIdx.x;
	int global_j = blockIdx.y * blockDim.y + threadIdx.y;

	for (int j = global_j; j < cuda::geometry.ysz; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::geometry.xsz - 1; i += blockDim.x * gridDim.x)
		{


			NUMERIC_TYPE h0 = U.H[j * cuda::geometry.xsz + i];
			NUMERIC_TYPE h1 = U.H[j * cuda::geometry.xsz + i + 1];
			NUMERIC_TYPE& qptr = U.Qx[j * (cuda::geometry.xsz + 1) + i + 1];

			qptr = C(0.0);

			if (h0 > cuda::solver_params.DepthThresh || h1 > cuda::solver_params.DepthThresh)
			{
				qptr = CalcFPQxAcc_(i, j, U, DEM, manning, friction);
			}
		}
	}
}


__global__ void
__launch_bounds__(CUDA_BLOCK_SIZE)
FloodplainQy
(

	Flow U,
	NUMERIC_TYPE* DEM,

	NUMERIC_TYPE* manning,

	bool friction
	)
{

	int global_i = blockIdx.x * blockDim.x + threadIdx.x;
	int global_j = blockIdx.y * blockDim.y + threadIdx.y;

	for (int j = global_j; j < cuda::geometry.ysz - 1; j += blockDim.y * gridDim.y)
	{
		for (int i = global_i; i < cuda::geometry.xsz; i += blockDim.x * gridDim.x)
		{


			NUMERIC_TYPE h0 = U.H[j * cuda::geometry.xsz + i];
			NUMERIC_TYPE h1 = U.H[(j + 1) * cuda::geometry.xsz + i /*+ 1*/];
			NUMERIC_TYPE& qptr = U.Qy[(j + 1) * (cuda::geometry.xsz + 1) + i /*+ 1*/];

			qptr = C(0.0);

			if (h0 > cuda::solver_params.DepthThresh || h1 > cuda::solver_params.DepthThresh)
			{
				qptr = CalcFPQyAcc_(i, j, U, DEM, manning, friction);
			}
		}
	}
}


}
}
}

lis::cuda::acc::Solver::Solver
(
	Flow& U,
	NUMERIC_TYPE* DEM,

	NUMERIC_TYPE* manning,
	Geometry& geometry,
	PhysicalParams& physical_params,
	dim3 grid_size
)
:
Uold(U1),
U(U1),
DEM(DEM),
manning(manning),
grid_size(grid_size),
maxH(geometry)
{
	Flow::allocate_device(U1, geometry);
//    Flow::allocate_device(U1, geometry);


	Flow::copy(U1, U, geometry);

	friction = (physical_params.manning > C(0.0)) || (manning != nullptr);
}

void lis::cuda::acc::Solver::zero_ghost_cells()
{
	zero_ghost_cells_north_south<<<1, CUDA_BLOCK_SIZE>>>(U);
	zero_ghost_cells_east_west<<<1, CUDA_BLOCK_SIZE>>>(U);
}

void lis::cuda::acc::Solver::update_ghost_cells()
{
	lis::cuda::acc::update_ghost_cells<<<64, CUDA_BLOCK_SIZE>>>(
			Uold, DEM, manning);
}


void lis::cuda::acc::Solver::update_uniform_rain(NUMERIC_TYPE rain_rate)
{
	update_uniform_rain_func << <grid_size, cuda::block_size >> > (
		Uold, DEM, rain_rate);
}


void lis::cuda::acc::Solver::updateMaxFieldACC
(
	NUMERIC_TYPE t
)
{
	update_max_field_ACC << <grid_size, cuda::block_size >> > (Uold, t);
}



lis::cuda::acc::Flow& lis::cuda::acc::Solver::update_flow_variables
(
	MassStats* mass_stats
)
{


	update_flow_variables_x<<<grid_size, cuda::block_size>>>
			(Uold);

	update_flow_variables_y<<<grid_size, cuda::block_size>>>
			(Uold);



	return Uold;
}







lis::cuda::acc::Flow& lis::cuda::acc::Solver::d_U()
{
	return Uold;
}

void lis::cuda::acc::Solver::update_dt_per_element 
(
	NUMERIC_TYPE* dt_field
) const
{
	lis::cuda::acc::update_dt_per_element<<<grid_size, cuda::block_size>>>(
			dt_field, Uold);		
}





























void lis::cuda::acc::Solver::FloodplainQ()
{
  // Calculate Qx
  FloodplainQx<<<grid_size, cuda::block_size>>>
			(Uold, DEM, manning, friction);

  // Calculate Qy
  FloodplainQy<<<grid_size, cuda::block_size>>>
  			(Uold, DEM, manning, friction);

  return;	
}

void lis::cuda::acc::Solver::zero_thin_depth_slopes()
{

  return;	
}

lis::cuda::acc::Solver::~Solver()
{
	Flow::free_device(U1);
//	Flow::free_device(U2);

}
