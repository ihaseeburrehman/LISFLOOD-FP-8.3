#include "cuda_dem.cuh"
#include "cuda_sample.cuh"
#include "cuda_simulate.cuh"
#include "cuda_solver.cuh"
#include "cuda_util.cuh"
#include "cuda_dg2_dem.cuh"
#include "cuda_dg2_flow.cuh"
#include "cuda_dg2_simulate.cuh"
#include "cuda_dg2_snapshot.cuh"
#include "cuda_dg2_solver.cuh"
#include "cuda_max_field.cuh"
#include "cuda_rain.cuh"
#include "cuda_stats.cuh"
#include "cuda_unifiedallocator.cuh"
#include "../rain/rain.h"
#include "io.h"
#include <cstdio>
#include <ctime>
#include "../../rain/Rainfall_uniform.h"


void lis::cuda::dg2::Simulation::run
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	::Solver& solver,
	int verbose
)
{
	if (verbose) print_device_info();
	dim3 grid_size(128, 128);

	Geometry geometry;
	int pitch = 0;
	int offset = 0;

	Topography DEM(filenames.demfilename, geometry, pitch, offset,
			pars.nodata_elevation, states.acceleration, verbose); 

	update_geometry(pars, geometry);
	cuda::copy_to_symbol(cuda::geometry, &geometry, sizeof(Geometry));
	cuda::copy_to_symbol(cuda::pitch, &pitch, sizeof(int));

	DeviceTopography d_DEM;
	d_DEM.initialise(DEM, geometry);

	cuda::sync();

    ::DynamicRain<cuda::UnifiedAllocator<NUMERIC_TYPE>> rain(
			filenames.dynamicrainfilename, verbose);
	cuda::copy_to_symbol(cuda::rain_geometry, &rain.geometry(), sizeof(Geometry));

	Flow U;
	Flow::allocate_pinned(U, geometry);
	initialise_H(U.H, U.H1x, U.H1y, filenames.startfilename, states, DEM._0,
			geometry, pitch, offset, verbose);
	initialise_discharge(U.HU, U.HV, filenames.startfilename, states,
			geometry, pitch, offset, verbose);

	NUMERIC_TYPE* manning = nullptr;
	NUMERIC_TYPE* d_manning = nullptr;
	if (strlen(filenames.nfilename) > 0)
	{
		manning = lis::GhostRaster::allocate(geometry);
		initialise_manning(manning, filenames.nfilename, geometry,
				pitch, offset, verbose);
		d_manning = lis::cuda::GhostRaster::allocate_device(geometry);
		GhostRaster::copy(d_manning, manning, geometry);
	}

	PhysicalParams physical_params(pars, solver);
	cuda::copy_to_symbol(cuda::physical_params, &physical_params,
			sizeof(PhysicalParams));

	SolverParams solver_params(pars, solver, states); 
	cuda::copy_to_symbol(cuda::solver_params, &solver_params,
			sizeof(SolverParams));

	Rainfall_uniform   rainfall_uniform(C(0.0));

	BoundCs boundCs;
	load_boundaries(filenames, states, pars, boundCs,
			verbose);
	cuda::Boundary::initialise(cuda::boundaries, boundCs, pitch, offset);

	auto solve = unique_ptr<cuda::Solver<Flow>>(new Solver(U, d_DEM, d_manning,
				geometry, physical_params, pars.limit_slopes, grid_size));

	DynamicTimestep<Flow> dynamic_dt(cuda::dt, geometry, solver.InitTstep,
			states.adaptive_ts, *solve);
	solver.MinTstep = solver.InitTstep;

	StatsCollector stats_collector(geometry, solver.DepthThresh, states.acceleration);
	MaxField max_field(filenames.resrootname, geometry, pitch, offset,
			grid_size, pars.output_precision);

	Stats stats(filenames.resrootname, pars.MassInt, pars.MassTotal, states.checkpoint, solver.t);
	stats.write_header(solver.t);

	Sampler sampler(cuda::sample_points, cuda::sample_buf_idx, verbose);
	if (states.save_stages == ON)
	{
		sampler.load_sample_points(filenames.stagefilename, geometry,
				pitch, offset);
		sampler.open_stage_file(filenames.resrootname, states.checkpoint, solver.t);
		sampler.write_stage_header(DEM._0, filenames.stagefilename, states.checkpoint, solver.t);

		if (states.voutput_stage == ON)
		{
			sampler.open_gauge_file(filenames.resrootname, states.checkpoint, solver.t);
			sampler.write_gauge_header(DEM._0, filenames.gaugefilename, solver.t);
		}
	}

	Snapshot snapshot(filenames.resrootname, pars.SaveInt,
			pars.SaveTotal, pars.SaveNo, DEM, U, geometry, pitch, offset,
			verbose, pars.output_precision);
	if (states.save_elev == ON)
	{
		snapshot.enable_elevation_writer();
	}
	if (states.voutput == ON)
	{
		snapshot.enable_velocity_writer(solver.DepthThresh);
	}
	if (states.save_Qs == ON)
	{
		snapshot.enable_discharge_writer();
	}

	Flow& d_U = solve->d_U();

	time_t loop_start;                                                          
	time(&loop_start);

	while (solver.t < solver.Sim_Time)
	{

		if (rain.enabled())
		{
			rain.update_time(solver.t);
			lis::cuda::DynamicRain::update<<<grid_size, CUDA_BLOCK_SIZE>>>(
					d_DEM._0, d_U.H, rain.data());
		}

		if (states.rainfall == ON) {

			rainfall_uniform.update_rainfall(filenames.rainfilename, solver.t);

			solve->update_uniform_rain(rainfall_uniform.value);
		}

		stats_collector.zero_instantaneous_mass();
		cuda::Boundary::update_time_series(solver.t);
		cuda::Boundary::update_point_sources(d_U.H, DEM._0,
				stats_collector.instantaneous_mass());
		stats_collector.accumulate_mass();
		stats_collector.zero_instantaneous_mass();
		solve->update_ghost_cells();
		dynamic_dt.update_dt();
		solver.MinTstep = FMIN(solver.MinTstep, cuda::dt);

		if (verbose) printf("t=%f\tdt=%f\n", solver.t, cuda::dt);
		solver.t += cuda::dt;
		solver.itCount++;

		d_U = solve->update_flow_variables(
				stats_collector.instantaneous_mass());
		if (pars.drain_nodata)
		{
			cuda::Boundary::drain_nodata_water(d_U.H, d_DEM._0, 
					stats_collector.instantaneous_mass(), grid_size);
		}
		stats_collector.accumulate_mass();
		max_field.update(d_U.H);

		if (stats.need_to_write(solver.t))
		{
			sampler.sample(d_U.H, d_U.HU, d_U.HV, solver.t);
			solve->zero_ghost_cells();
			cuda::sync();

			stats.write(stats_collector.create_entry(solver, d_U.H));
			sampler.write_if_buffer_full();
		}

		snapshot.write_if_needed(d_U, solver.t, states.call_gzip);
	}

	time_t loop_end;
	time(&loop_end);
	NUMERIC_TYPE seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);   

	sampler.write();
	snapshot.wait();
	max_field.write();

	d_DEM.free_device();
	delete[] manning;
	cuda::free_device(d_manning);
	Flow::free_pinned(U);
	cuda::Boundary::free_device(boundaries);
}

void lis::cuda::dg2::Simulation::initialise_H
(
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* H1x,
	NUMERIC_TYPE* H1y,
	const char* filename,
	States& states,
	NUMERIC_TYPE* DEM,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	if (states.startfile == ON)
	{
		cuda::Simulation::initialise_H(H, filename, states, DEM, geometry,
				pitch, offset, verbose);
		initialise_H_slope(H1x, filename, "1x", geometry, pitch, offset,
				verbose);
		initialise_H_slope(H1y, filename, "1y", geometry, pitch, offset,
				verbose);
	}
}

void lis::cuda::dg2::Simulation::initialise_H_slope
(
	NUMERIC_TYPE* array,
	const char* filename,
	const char* suffix,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	char h_slope[256]; 
	strcpy(h_slope, filename);
	strcat(h_slope, suffix); 

	StartFile::load(h_slope, array, geometry, pitch, offset, verbose);
}
