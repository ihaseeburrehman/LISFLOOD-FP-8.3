#pragma once
#include "cuda_dg2_dem.cuh"
#include "cuda_dg2_flow.cuh"
#include "cuda_dg2_slope_limit.cuh"
#include "cuda_solver.cuh"

namespace lis
{
namespace cuda
{
namespace dg2
{

extern __managed__ NUMERIC_TYPE maxH;

class Solver : public cuda::Solver<Flow>
{
public:
	Solver
	(
		Flow& U,
		DeviceTopography& DEM,
		NUMERIC_TYPE* manning,
		Geometry& geometry,
		PhysicalParams& physical_params,
		bool limit_slopes,
		dim3 grid_size
	);

	void FloodplainQ(); 

	void zero_ghost_cells();

	void update_ghost_cells();
	
	void update_uniform_rain(NUMERIC_TYPE  rain_rate);
	
	void zero_thin_depth_slopes();

	void updateMaxFieldACC(NUMERIC_TYPE t);

	Flow& update_flow_variables
	(
		MassStats* mass_stats
	);

	Flow& d_U();

	void update_dt_per_element
	(
		NUMERIC_TYPE* dt_field
	) const;

	~Solver();

private:
	void update_ghost_cells
	(
		Flow& U
	);
	
	void zero_thin_depth_slopes
	(
		Flow& U
	);	

	void zero_perimeter_slopes
	(
		DeviceTopography& DEM,
		Flow& flow
	);

	Flow U1;
	Flow Ux;
	Flow Uint;
	Flow& Uold;
	DeviceTopography& DEM;
	NUMERIC_TYPE* manning;
	Slopes slopes;
	MaxH maxH;
	bool friction;
	bool limit_slopes = false;
	const dim3 grid_size;
};

}
}
}
