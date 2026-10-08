#include "cuda_dem.cuh"
#include "cuda_sample.cuh"
#include "cuda_simulate.cuh"
#include "cuda_acc_flow.cuh"
#include "cuda_acc_simulate.cuh"
#include "cuda_acc_snapshot.cuh"
#include "cuda_acc_solver.cuh"
#include "cuda_rain.cuh"
#include "cuda_max_field.cuh"
#include "cuda_stats.cuh"
#include "cuda_unifiedallocator.cuh"
#include "../rain/rain.h"
#include "io.h"
#include <cstdio>
#include <ctime>
#include "stats.h"
#include "../../rain/Rainfall_uniform.h"

//#include "../acc_nugrid/CHECK_CUDA_ERROR.cuh" // use for debug

void lis::cuda::acc::Simulation::run
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

	NUMERIC_TYPE* DEM = Topography::load(filenames.demfilename, geometry,
			pitch, offset, pars.nodata_elevation, states.acceleration, verbose); 
	
	//CHECK_CUDA_ERROR(cuda::peek());
	//CHECK_CUDA_ERROR(cuda::sync());

	update_geometry(pars, geometry); 

	cuda::copy_to_symbol(cuda::geometry, &geometry, sizeof(Geometry)); 

	cuda::copy_to_symbol(cuda::pitch, &pitch, sizeof(int));

	NUMERIC_TYPE* d_DEM = GhostRaster::allocate_device_H(geometry); 

	GhostRaster::copy_H(d_DEM, DEM, geometry);

//	cuda::sync();

    ::DynamicRain<cuda::UnifiedAllocator<NUMERIC_TYPE>> rain( 
			filenames.dynamicrainfilename, verbose);

	cuda::copy_to_symbol(cuda::rain_geometry, &rain.geometry(), sizeof(Geometry));

	acc::Flow U;
	
	acc::Flow::allocate_pinned(U, geometry);
	initialise_H(U.H, filenames.startfilename, states, DEM, 
			geometry, pitch, offset, verbose);

	nullify_max(U.maxHtm, U.initHtm, U.totalHtm, geometry, pitch, offset, verbose);

	NUMERIC_TYPE* manning = nullptr;
	NUMERIC_TYPE* d_manning = nullptr;
	if (strlen(filenames.nfilename) > 0) 
	{
		manning = lis::GhostRaster::allocate_H(geometry);
		initialise_manning(manning, filenames.nfilename, geometry,
				pitch, offset, verbose);
		d_manning = lis::cuda::GhostRaster::allocate_device_H(geometry);
		GhostRaster::copy_H(d_manning, manning, geometry);
	}


	BoundCs boundCs;
	load_boundaries(filenames, states, pars, boundCs, verbose);

	Rainfall_uniform   rainfall_uniform(C(0.0));

//	MaxField max_field(filenames.resrootname, geometry, pitch, offset,
//		grid_size, states.acceleration, pars.output_precision);

	

	StatsCollector stats_collector(geometry, solver.DepthThresh, states.acceleration);


	//Load checkpointed data if this job has been restarted
	if (states.checkpoint == ON)
	{
		
		read_checkpoint(filenames, states, pars, solver,
			U.H, U.Qx, U.Qy, U.maxH, U.totalHtm, U.maxHtm, U.initHtm, stats_collector, verbose);

		if (verbose == ON) printf(" - checkpoint output file: %s\n", filenames.checkpointfilename);
	}

	PhysicalParams physical_params(pars, solver); 

	cuda::copy_to_symbol(cuda::physical_params, &physical_params,
			sizeof(PhysicalParams));

	SolverParams solver_params(pars, solver, states); 

	cuda::copy_to_symbol(cuda::solver_params, &solver_params,
			sizeof(SolverParams));

	cuda::Boundary::initialise(cuda::boundaries, boundCs, pitch, offset); 

	auto solve = unique_ptr<cuda::Solver<Flow>>(new Solver( 
				U, d_DEM, d_manning, geometry,
				physical_params, grid_size));				

	DynamicTimestepACC<acc::Flow> dynamic_dt(cuda::dt, geometry, solver.InitTstep, 
			states.adaptive_ts, *solve); 


	//start simulation
	time(&solver.time_start);

	solver.itrn_time_now = solver.itrn_time;
    solver.MinTstep = solver.InitTstep;

//	StatsCollector stats_collector(geometry, solver.DepthThresh);

	Stats stats(filenames.resrootname, pars.MassInt, pars.MassTotal, states.checkpoint, solver.t);

	stats.write_header(solver.t);

	Sampler sampler(cuda::sample_points, cuda::sample_buf_idx, verbose);

	if (states.save_stages == ON)
	{
		sampler.load_sample_points(filenames.stagefilename, geometry,
				pitch, offset);

		sampler.open_stage_file(filenames.resrootname, states.checkpoint, solver.t);

		sampler.write_stage_header(DEM, filenames.stagefilename, states.checkpoint, solver.t);

		if (states.voutput_stage == ON)
		{
			sampler.open_gauge_file(filenames.resrootname, states.checkpoint, solver.t);

			sampler.write_gauge_header(DEM, filenames.gaugefilename, solver.t);
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
			lis::cuda::DynamicRain::updateACC << <grid_size, CUDA_BLOCK_SIZE >> > (
				d_DEM, d_U.H, rain.data());
		}

		stats_collector.zero_instantaneous_mass(); 

		dynamic_dt.update_dt_ACC(); 

		solve->FloodplainQ();

		cuda::Boundary::update_time_series(solver.t); 

        solve->update_ghost_cells(); 

		if (pars.drain_nodata)
		{

			cuda::Boundary::drain_nodata_waterACC(d_U.H, d_DEM,
					stats_collector.instantaneous_mass(), grid_size); 
		}

		if (states.rainfall == ON) {

			rainfall_uniform.update_rainfall(filenames.rainfilename, solver.t);

			solve->update_uniform_rain(rainfall_uniform.value);
		}

		stats_collector.accumulate_mass();

		stats_collector.zero_instantaneous_mass();

		cuda::Boundary::update_H(d_U.H, d_DEM, d_U.Qx, d_U.Qy, 
				stats_collector.instantaneous_mass(), grid_size, states.drychecking);

		stats_collector.accumulate_mass();
	

		d_U = solve->update_flow_variables( 
				stats_collector.instantaneous_mass());
		
		cuda::sync();
	
		solver.MinTstep = FMIN(solver.MinTstep, cuda::dt);

		if (verbose) printf("t=%f\tdt=%f\n", solver.t, cuda::dt);
		solver.t += cuda::dt; 
		solver.itCount++;

		//if ((Statesptr->reset_timeinit == ON) & (Solverptr->t > Parptr->reset_timeinit_time))
		//{ //reset the time of initial inundation if called for in parameter file
		//	Statesptr->reset_timeinit = OFF;

		//	//#pragma omp parallel for private(j,ptr) // irrelevent benefit JCN
		//	for (j = 0; j < Parptr->ysz; j++) for (i = 0; i < Parptr->xsz; i++)
		//	{
		//		ptr = i + j * Parptr->xsz;
		//		Arrptr->initHtm[ptr] = (NULLVAL);
		//	}
		//	if (verbose == ON) printf("\n Time of initial inundation reset \n");
		//}


//		max_field.updateACC(d_U.H);
		
		// Update maxH, maxHtm, totalHtm and initHtm at the mass interval if mint_hk is specifed in the .par file OR at every time step if not
		if (stats.need_to_write(solver.t) || states.mint_hk == OFF)
		{
			solve->updateMaxFieldACC(solver.t);
		}

		snapshot.write_if_needed(d_U, solver.t, states.call_gzip);
		pars.SaveTotal = snapshot.next_save;
		pars.SaveNo = snapshot.counter;

		if (stats.need_to_write(solver.t))
		{
//			solve->updateMaxFieldACC(solver.t);

			sampler.sample_ACC(d_U.H, d_U.Vx, d_U.Vy, solver.t); 

//			lis::StatsEntry Stats_Entry;

			lis::StatsEntry Stats_Entry = stats_collector.create_entry(solver, d_U.H);

			stats.write(Stats_Entry);

			pars.MassTotal = stats.next_save;

			sampler.write_if_buffer_full();

			// Checkpointing
			if (states.checkpoint == ON)
			{
				//iteration time
				time(&solver.time_check);
				solver.itrn_time_now = solver.itrn_time + (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start);

				Flow::copy(U, d_U, geometry);

				if (solver.itrn_time_now >= pars.nextcheck)
				{
//					NUMERIC_TYPE SaveTotal, MassTotal;
//					int SaveNo;

					
					write_checkpoint(filenames, states, pars, solver,
						U.H, U.Qx, U.Qy, U.maxH, U.totalHtm, U.maxHtm, U.initHtm, Stats_Entry, cuda::dt, verbose);

					pars.nextcheck = solver.itrn_time_now + (pars.checkfreq * 3600);
				}
			}
		}



		



		// If requested on command line, check whether we should kill this simulation...
		if (states.killsim == ON)
		{
			//iteration time
			time(&solver.time_check);
			solver.itrn_time_now = solver.itrn_time + (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start);
			// check if we have reached the kill time
			if (solver.itrn_time_now >= pars.killsim_time) {
				if (verbose == ON) printf("Simulation kill time reached... ");
				break;
			}
		}

	}

	time_t loop_end;
	time(&loop_end);
	NUMERIC_TYPE seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);  

//	sampler.write();
//	snapshot.wait();
//	max_field.write();

	//output max files
	snapshot.write_max_files(d_U);




	// Must be last because maxH changed to Maximum elevation !!!
	// Write maximum elevation
	//for (i = 0; i < pars.xsz; i++) for (j = 0; j < pars.ysz; j++)
	//{
	//	ptr = i + j * pars.xsz;
	//	if (d_U.maxH[ptr] > 1e-3) d_U.maxH[ptr] += d_DEM[ptr]; else d_U.maxH[ptr] = NULLVAL;
	//}
 //   

	//write_ascfile(filenames.resrootname, -1, ".mxe", d_U.maxH, d_DEM, 0, states, pars);

	//if (verbose == ON) printf("Finished.\n\n");



	delete[] DEM;
	delete[] manning;
	cuda::free_device(d_DEM);
	cuda::free_device(d_manning);
	Flow::free_pinned(U);
	cuda::Boundary::free_device(boundaries);
}
