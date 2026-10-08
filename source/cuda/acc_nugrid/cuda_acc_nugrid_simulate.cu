// RC's suggestion to manage Intellisense
#ifdef __INTELLISENSE__
#ifndef __CUDACC__
#define __CUDACC__
#endif
#endif

#include "cuda_acc_nugrid_simulate.cuh"
#include "io.h"
#include <cstdio>
#include <ctime>

//-----------------Custom headers-----------------//

// Kernels
#include "generate_all_morton_codes.cuh"
#include "copy_finest_coefficients.cuh"
//#include "zero_details.cuh"
#include "count_neighbours.cuh"
#include "find_nonuniform_neighbours.cuh"
#include "get_compaction_flags.cuh"
#include "count_interfaces_per_neighbours.cuh"
#include "sort_nghbr_assem_sol_row_major.cuh"
#include "init_neighbours.cuh"
#include "init_interfaces.cuh"
#include "find_interfaces.cuh"
#include "load_interface_q_vol.cuh"
#include "load_neighbour_h.cuh"
#include "compute_q.cuh"
#include "update_h.cuh"
#include "calculate_dt.cuh"
#include "init_h.cuh"

// Kernel wrappers
#include "get_modal_values.cuh"
#include "sort_finest_scale_coefficients_z_order.cuh"
#include "get_max_scale_coefficients.cuh"
#include "preflag_topo.cuh"
#include "get_reg_tree.cuh"
#include "reverse_z_order_assembled_solution.cuh"
#include "reverse_z_order_act_idcs.cuh"
#include "compaction.cuh"
#include "get_dt_CFL.cuh"
#include "read_all_bdy_conds.cuh"

// Input/output
#include "write_soln_vtk.cuh"
#include "write_gauge_point_data.cuh"
#include "read_gauge_points.cuh"
#include "read_point_srcs.cuh"
#include "read_num_stage_points.h"
//#include "read_num_all_bdy_cells.h"
#include "read_num_point_srcs.h"
#include "init_q.cuh"

// Helper functions
#include "get_lvl_idx.cuh"
#include "preflag_details.cuh"
#include "project_assem_sol.cuh"
#include "../../rain/rain.h"
#include "rain.cuh"
#include "unifiedallocator.cuh"
#include "write_mass_data.cuh"
#include "RainfallUniform.h"
#include "rain_uniform.cuh"
#include "update_max.cuh"
#include "write_max_maps.cuh"
#include "init_max.cuh"
#include "write_chkpnt.cuh"
#include "read_chkpnt.cuh"

void lis::cuda::acc_nugrid::Simulation::run
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	::Solver& solver,
	int verbose
)
{
	if (verbose) print_device_info();
	dim3 grid_size(128, 128); //TOREMOVE

//	Geometry geometry;

    NUMERIC_TYPE Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;

	// =========================================================== //
	// INITIALISATION OF VARIABLES AND INSTANTIATION OF STRUCTURES //
	// =========================================================== //

	// Variables
	int mesh_dim = 1 << solver.L;
	int interface_dim = mesh_dim + 1;

	int num_finest_elems      = mesh_dim * mesh_dim;
	int num_blocks_finest     = get_num_blocks(num_finest_elems, THREADS_PER_BLOCK_MRA); 
	int num_threads_traversal = num_finest_elems / 4;
	int num_blocks_traversal  = get_num_blocks(num_threads_traversal, THREADS_PER_BLOCK_MRA);
	int num_all_elems         = get_lvl_idx(solver.L + 1);
	int num_details           = get_lvl_idx(solver.L);
	int num_blocks_details    = get_num_blocks(num_details, THREADS_PER_BLOCK_MRA);
	int num_blocks_sol        = 0;
	int num_blocks_itfc       = 0;
	int num_blocks_nghbrs     = 0;
	int num_blocks_all        = get_num_blocks(num_all_elems, THREADS_PER_BLOCK_MRA);

	index_1D finest_lvl_idx   = get_lvl_idx(solver.L);
	clock_t end               = clock(); //TOREMOVE
	bool    non_uniform_n     = false;

	if (strlen(filenames.nfilename) > 0) non_uniform_n = true;

//	CHECK_CUDA_ERROR(peek());
//	CHECK_CUDA_ERROR(sync());

	BoundCs boundCs;

	LoadBCs(&filenames, &states, &pars, &boundCs, verbose);
	LoadBCVar(&filenames, &states, &pars, &boundCs, nullptr, nullptr, nullptr, verbose);

	// Structures
	Maxes           maxes         = { C(0.0) }; 
	GaugePoints     gauge_points  ( read_num_stage_points(filenames.stagefilename));
//	Boundaries      boundaries    ( read_num_all_bdy_cells(filenames.bcifilename, pars));
    Boundaries      boundaries(filenames, pars);
	PointSources    point_sources ( read_num_point_srcs(filenames.bcifilename));
	RainfallUniform rainfall_uniform (C(0.0));

	AssembledSolution d_assem_sol(num_finest_elems, non_uniform_n);
	AssembledSolution d_buf_assem_sol(num_finest_elems, non_uniform_n);
	AssembledSolution d_plot_assem_sol(num_finest_elems, non_uniform_n);
	ScaleCoefficients d_scale_coeffs(num_all_elems, non_uniform_n);
	Details           d_details(num_details, non_uniform_n, states.startfile);
	CompactionFlags   d_compaction_flags(num_finest_elems);
	AssembledSolution d_nghbr_assem_sol(num_finest_elems, non_uniform_n);
	
	read_gauge_points(filenames.stagefilename, pars, gauge_points);
//	read_all_bdy_conds(filenames.bcifilename, boundaries, pars);
	read_point_srcs(filenames.bcifilename, filenames.bdyfilename, pars, solver.t, point_sources, boundCs);

	// Bytesizes
	size_t bytes_morton = num_finest_elems * sizeof(MortonCode);
//	size_t bytes_details = num_details * sizeof(NUMERIC_TYPE);
	size_t bytes_soln = num_finest_elems * sizeof(NUMERIC_TYPE);
//	size_t bytes_sig_details = num_details * sizeof(bool); 

	// Arrays
	MortonCode* d_morton_codes = (MortonCode*)malloc_device(bytes_morton);
	MortonCode* d_sorted_morton_codes = (MortonCode*)malloc_device(bytes_morton);
	MortonCode* d_indices = (MortonCode*)malloc_device(bytes_morton);
	MortonCode* d_reverse_z_order = (MortonCode*)malloc_device(bytes_morton);
	NUMERIC_TYPE* d_eta_temp = (NUMERIC_TYPE*)malloc_device(bytes_soln); //TOREMOVE
	bool* d_sig_details = preflag_details(boundaries, point_sources, gauge_points, num_details, solver.L); 
	NUMERIC_TYPE* d_dt_CFL = (NUMERIC_TYPE*)malloc_device(bytes_soln);
	NUMERIC_TYPE* maxH = (NUMERIC_TYPE*)malloc_device(bytes_soln);
	NUMERIC_TYPE* totalHtm = (NUMERIC_TYPE*)malloc_device(bytes_soln);
	NUMERIC_TYPE* maxHtm = (NUMERIC_TYPE*)malloc_device(bytes_soln);
	NUMERIC_TYPE* initHtm = (NUMERIC_TYPE*)malloc_device(bytes_soln);

	// =========================================================== //

	::DynamicRain<UnifiedAllocator<NUMERIC_TYPE>> rain(
		filenames.dynamicrainfilename, verbose);

	// ================================ //
	// PREPROCESSING BEFORE SOLVER LOOP //
	// ================================ //

	init_max << <num_blocks_finest, THREADS_PER_BLOCK_MRA >> >
	(
		d_plot_assem_sol,
		maxH,
		totalHtm,
		maxHtm,
		initHtm
	);

	get_modal_values
	(
		d_buf_assem_sol,
		mesh_dim,
		interface_dim,
		filenames,
        states,
		pars.nodata_elevation
	);

	generate_all_morton_codes<<<num_blocks_finest, THREADS_PER_BLOCK_MRA>>>
	(
		d_morton_codes,
		d_indices,
		mesh_dim
	);

	sort_finest_scale_coefficients_z_order
	(
		d_morton_codes, // unsorted
		d_sorted_morton_codes, // sorted
		d_buf_assem_sol, // unsorted // row majored
		d_assem_sol, // sorted // zorder
		d_indices, // unsorted
		d_reverse_z_order, // sorted
		non_uniform_n
	);

	copy_finest_coefficients<<<num_blocks_finest, THREADS_PER_BLOCK_MRA>>>
	(
		d_assem_sol, // sorted // has h 
		d_scale_coeffs, // sorted // has eta
		finest_lvl_idx,
		non_uniform_n
	);

	maxes = get_max_scale_coefficients(d_assem_sol); 

	preflag_topo /// flag cells with significant topography detail
	(
		d_scale_coeffs, /// in /// updated
		d_details, /// out /// updated only topo
		d_sig_details, /// in /// updated
		maxes,
		solver.epsilon,
		solver.L,
		non_uniform_n,
        states.startfile
	);

	// ================================ //

    get_reg_tree
    (
    	d_sig_details,
    	solver.L
    );

    if (non_uniform_n) {
    	traverse_tree_of_sig_details_with_n << <num_blocks_traversal, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_sig_details,
    		d_scale_coeffs,
    		d_buf_assem_sol, // out // changes it from row major to non-compacted z order //only the assembled grid not the whole hierarchy
    		num_threads_traversal,
    		solver.L
    	);
    } else {
    
    	traverse_tree_of_sig_details << <num_blocks_traversal, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_sig_details,
    		d_scale_coeffs,
    		d_buf_assem_sol, // out // non-compacted z order //only the assembled grid not the whole hierarchy
    		num_threads_traversal,
    		solver.L
    	);
    }

    reverse_z_order_act_idcs
    (
    	d_reverse_z_order, // d_keys_in
    	d_indices,         // d_keys_out
    	d_buf_assem_sol,   // d_values_in // non-compacted z order //only the assembled grid not the whole hierarchy
    	d_assem_sol,       // d_values_out // non-compacted row major // to be copied into another array to use for find_nonuniform_neighbours
    	num_finest_elems
    );

    get_compaction_flags << <num_blocks_finest, THREADS_PER_BLOCK_MRA >> >
    (
    	d_buf_assem_sol,
    	d_compaction_flags,
    	num_finest_elems
    );

    compaction
    (
    	d_assem_sol, // out // compacted z order
    	d_buf_assem_sol, // in // non-compacted z order
    	d_compaction_flags,
    	num_finest_elems,
    	non_uniform_n,
        states.startfile
    );

    // GRID DIMENSIONS BASED ON ASSEMBLED SOLUTION LENGTH //
    num_blocks_sol = get_num_blocks(d_assem_sol.length, THREADS_PER_BLOCK_MRA);
    
    if (!states.startfile) {
        init_h << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> >
            (
                d_assem_sol
            );
    }

    sort_nghbr_assem_sol_row_major // exactly the same as reverse_z_order_act_idcs
    (
    	d_reverse_z_order,
    	d_indices,
    	d_buf_assem_sol, // in // non-compacted z order
    	d_nghbr_assem_sol, // out // non-compacted row major 
    	num_finest_elems
    );

    count_neighbours << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> >
    (
    	d_assem_sol, // compacted z order
    	d_nghbr_assem_sol,
    	pars,
    	solver,
    	boundaries
    );

    NonUniformNeighbours d_non_uniform_nghbrs = init_neighbours(d_assem_sol, non_uniform_n);

    find_nonuniform_neighbours << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> >
    (
    	d_assem_sol, // compacted z order
    	d_nghbr_assem_sol, // non-compacted row major
    	pars,
    	solver,
    	d_non_uniform_nghbrs,
    	boundaries
    );

    num_blocks_nghbrs = get_num_blocks(d_non_uniform_nghbrs.length, THREADS_PER_BLOCK_MRA);
    
    init_q << <num_blocks_nghbrs, THREADS_PER_BLOCK_MRA >> >
    (
    	d_non_uniform_nghbrs
    );


    //Load checkpointed data if this job has been restarted
    if (states.checkpoint == ON)
    {
        read_chkpnt(d_assem_sol, d_non_uniform_nghbrs, filenames, states, pars, solver,
            maxH, totalHtm, maxHtm, initHtm, non_uniform_n, num_finest_elems, verbose);

        if (verbose == ON) printf(" - checkpoint output file: %s\n", filenames.checkpointfilename);
    }

    //mass balance

    char fullpath[255];

    sprintf(fullpath, "%s%s", filenames.resrootname, ".mass");

    FILE* fp;

    if (states.checkpoint == ON && solver.t > 0) { //if this is a checkpointed job, we only need to amend the .mass file
        fp = fopen(fullpath, "a");
    }
    else {
        fp = fopen(fullpath, "w");
    }
    if (fp != NULL)
    {
        if (states.checkpoint == ON && solver.t > 0)
        {
            // make a note in the mass file that this is a restart point - user can then edit the overlap out if they want a continuous mass file record.
            fprintf(fp, "####################################################### Checkpoint restart ########################################################\n");
            fprintf(fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-(Inf+Evap)\n");
            fflush(fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
        }
    }
    else
    {
        if (verbose == ON)
        {
            printf("Unable to open mass balance file: %s", fullpath);
            exit(0);
        }
    }

    fclose(fp);


    write_mass_data
    (
        filenames.resrootname,
        gauge_points,
        mesh_dim,
        solver,
        states,
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0),
        C(0.0)
    );


    //stage output file
    if (states.save_stages == ON) {

        sprintf(fullpath, "%s%s", filenames.resrootname, ".stage");

        if (states.checkpoint == ON && solver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
            fp = fopen(fullpath, "a");
        }
        else {
            fp = fopen(fullpath, "w");
        }

        if (fp != NULL)
        {
            if (states.checkpoint == ON && solver.t > 0)
            {
                fprintf(fp, "####################################################### Checkpoint restart ########################################################\n");
                fflush(fp);
            }
        }
        else
        {
            if (verbose == ON) printf("Unable to open stage output file: %s", fp);
            states.save_stages = OFF;
        }
    }

    fclose(fp);

    //velocity output file
    if (states.save_stages == ON && states.voutput_stage == ON) {

        sprintf(fullpath, "%s%s", filenames.resrootname, ".velocity");

        if (states.checkpoint == ON && solver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
            fp = fopen(fullpath, "a");
        }
        else {
            fp = fopen(fullpath, "w");
        }

        if (fp != NULL)
        {
            if (solver.t > 0)
            {
                fprintf(fp, "####################################################### Checkpoint restart ########################################################\n");
                fflush(fp);
            }
        }
        else
        {
            if (verbose == ON) printf("Unable to open stage output file: %s", fp);
            states.save_stages = OFF;
        }
    }

    fclose(fp);


    if (states.save_stages == ON)
    {
        d_plot_assem_sol = project_assem_sol
        (
            mesh_dim,
            d_sig_details,
            d_scale_coeffs,
            d_buf_assem_sol,
            solver.L,
            d_reverse_z_order,
            d_indices,
            d_assem_sol,
            d_plot_assem_sol
        );

        write_gauge_point_data
        (
            filenames.resrootname,
            d_plot_assem_sol,
            gauge_points,
            mesh_dim,
            filenames.stagefilename,
            pars,
            solver,
            states
        );

        if (states.voutput_stage == ON) {

            write_velocity_point_data
            (
                filenames.resrootname,
                d_plot_assem_sol,
                gauge_points,
                mesh_dim,
                filenames.stagefilename,
                pars,
                solver,
                states
            );
        }
    }


    count_interfaces_per_neighbours << <num_blocks_nghbrs, THREADS_PER_BLOCK_MRA >> >
    (
    	d_non_uniform_nghbrs,
    	d_assem_sol // compacted z order
    );

    NonUniformInterfaces d_non_uniform_itrfaces = init_interfaces(d_non_uniform_nghbrs);

    find_interfaces << <num_blocks_nghbrs, THREADS_PER_BLOCK_MRA >> >
    (
    	d_non_uniform_nghbrs,
    	d_assem_sol, // compacted z order
    	d_non_uniform_itrfaces
    );

    num_blocks_itfc = get_num_blocks(d_non_uniform_itrfaces.length, THREADS_PER_BLOCK_MRA);

    //start simulation
    time(&solver.time_start);

    if (verbose == ON)
    {
    	printf("\nStarting time steps: ");
    	fflush(stdout);
    }

    // Populating Tstep variables prior to start of simulation
    solver.itrn_time_now = solver.itrn_time;

    if (solver.t == 0)
	{
         solver.Tstep = solver.InitTstep;
         solver.MinTstep = solver.InitTstep;
	}
    if (verbose == ON) printf("acceleration mode\n\n");
    fflush(stdout);

    time_t loop_start;
    time(&loop_start);

	// main iteration loop
    while (solver.t < solver.Sim_Time)
    {
        if (pars.drain_nodata)
        {
            drain_nodata_water << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> > (
                d_assem_sol, pars, solver);
        }

        sync();

    	if (rain.enabled())
    	{
    		rain.update_time(solver.t);
    		update_rain << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> > (
    			d_assem_sol, rain.data(), rain.geometry().xsz, rain.geometry().tly, 
                rain.geometry().dx, rain.geometry().dy, pars, solver);
        }
    
    	if (states.rainfall == ON) {
    
    		rainfall_uniform.update_rainfall(filenames.rainfilename, solver.t);
    
    		update_uniform_rain << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> > (
    			d_assem_sol, pars, solver.Tstep, rainfall_uniform.value);
    	}
    
    	// Time step is reset to initial time step at the start of each FloodplainQ calculation
    	if (solver.t > C(0.0)) solver.Tstep = solver.InitTstep;

    	calculate_dt << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_assem_sol,
    		pars,
    		solver,
    		d_dt_CFL
    	);

        solver.Tstep = get_dt_CFL(d_dt_CFL, d_assem_sol.length);
    
    	point_sources.update_all_sources(filenames.bdyfilename, solver.t, boundCs);
    
    	boundaries.update_all_inlets(filenames.bdyfilename, solver.t/*, boundCs*/);

    	load_interface_q_vol << <num_blocks_itfc, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_non_uniform_nghbrs,
    		d_non_uniform_itrfaces
    	);

    	load_neighbour_h << <num_blocks_nghbrs, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_assem_sol,
    		d_non_uniform_nghbrs,
    		non_uniform_n
    	);

    	compute_q << <num_blocks_nghbrs, THREADS_PER_BLOCK_MRA >> >
    	(
            d_assem_sol,
            d_non_uniform_nghbrs,
            d_non_uniform_itrfaces,
            pars,
            solver,
            boundaries,
            non_uniform_n
    	);

    	update_h << <num_blocks_sol, THREADS_PER_BLOCK_MRA >> >
    	(
    		d_assem_sol,
    		d_non_uniform_nghbrs,
    		d_non_uniform_itrfaces,
    		pars,
    		solver,
    		boundaries,
    		point_sources
    	);

    	// Update t with final Tstep calculated
    	if (solver.t > C(0.0)) solver.MinTstep = getmin(solver.MinTstep, solver.Tstep);
    	solver.t += solver.Tstep;
    	solver.itCount += 1;
    
    	if (solver.t >= pars.MassTotal)
    	{
    		write_mass_data
    		(
    			filenames.resrootname,
    			gauge_points,
    			mesh_dim,
    			solver,
                states,
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0),
    			C(0.0)
    		);

            pars.MassTotal += pars.MassInt;
    
    		d_plot_assem_sol = project_assem_sol
    		(
    			mesh_dim,
    			d_sig_details,
    			d_scale_coeffs,
    			d_buf_assem_sol,
    			solver.L,
    			d_reverse_z_order,
    			d_indices,
    			d_assem_sol,
    			d_plot_assem_sol
    		);
    
    		update_max << <num_blocks_finest, THREADS_PER_BLOCK_MRA >> >
    		(
    			d_plot_assem_sol,
    			solver,
    			maxH,
    			totalHtm,
    			maxHtm,
    			initHtm
    		);
    
    		printf("Writing outputs at time %lf\n", solver.t);
    
            if (states.save_stages == ON) {
                write_gauge_point_data
                (
                    filenames.resrootname,
                    d_plot_assem_sol,
                    gauge_points,
                    mesh_dim,
                    filenames.stagefilename,
                    pars,
                    solver,
                    states
                );

                if (states.voutput_stage == ON) {
                    write_velocity_point_data
                    (
                        filenames.resrootname,
                        d_plot_assem_sol,
                        gauge_points,
                        mesh_dim,
                        filenames.stagefilename,
                        pars,
                        solver,
                        states
                    );
                }
            }

    		// Checkpointing
    		if (states.checkpoint == ON)
    		{
    			//iteration time
    			time(&solver.time_check);
    			solver.itrn_time_now = solver.itrn_time + (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start);
    
    			if (solver.itrn_time_now >= pars.nextcheck)
    			{
    				write_chkpnt(d_assem_sol, d_non_uniform_nghbrs, filenames, states, pars, solver,
    					maxH, totalHtm, maxHtm, initHtm, non_uniform_n, num_finest_elems, verbose);
    
    				pars.nextcheck = solver.itrn_time_now + (pars.checkfreq * 3600);
    			}
    		}
    	}
    
    	if (solver.t >= pars.SaveTotal) {
    
            time(&solver.time_check);
            Comp_time = (NUMERIC_TYPE)difftime(solver.time_check, solver.time_start) / 60;
            if (Comp_time != 0 && states.comp_out == ON) // only of t is not zero (can't divide by zero)
            {
                Model_Comp_Ratio = ((solver.t / 60) / Comp_time);
                Model_time_left = (solver.Sim_Time - solver.t) / 60;
                Est_Time_Fin = (Model_time_left / Model_Comp_Ratio);
                Est_Time_Tot = Comp_time + Est_Time_Fin;
                printf("T(mins): M: %.1" NUM_FMT", C: %.1" NUM_FMT", M/C: %.2" NUM_FMT", ETot: %.1" NUM_FMT", EFin: %.1" NUM_FMT"\n", (solver.t / C(60.0)), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
            }

//    		d_plot_assem_sol = project_assem_sol
//    		(
//    			mesh_dim,
//    			d_sig_details,
//    			d_scale_coeffs,
//    			d_buf_assem_sol,
//    			solver.L,
//    			d_reverse_z_order,
//    			d_indices,
//    			d_assem_sol,
//    			d_plot_assem_sol
//    		);
    
//    		write_raster_maps
//    		(
//    			filenames.resrootname,
//    			d_plot_assem_sol,
//    			mesh_dim,
//    			pars,
//                states.call_gzip,
//    			pars.output_precision
//    		);
    
            if (states.save_vtk == ON) {
                if (non_uniform_n) {
                    write_soln_vtk_with_n
                    (
                        filenames.resrootname,
                        d_assem_sol,
                        pars,
                        solver.L,
                        solver.DepthThresh,
                        states.call_gzip,
                        pars.output_precision
                    );
                }
                else {
                    write_soln_vtk
                    (
                        filenames.resrootname,
                        d_assem_sol,
                        pars,
                        solver.L,
                        solver.DepthThresh,
                        states.call_gzip,
                        pars.output_precision
                    );
                }
            }
            // update interval counter
            pars.SaveTotal += pars.SaveInt;
            pars.SaveNo += 1;
    	}
    }
    //END main ITERATIONS
    
    time_t loop_end;
    time(&loop_end);
    
    NUMERIC_TYPE seconds = difftime(loop_end, loop_start);
    printf("loop time %lf\n", seconds);
    
    //output max files
    write_max_maps
    (
    	filenames.resrootname,
    	mesh_dim,
    	pars,
    	maxH,
    	totalHtm,
    	maxHtm,
    	initHtm,
    	pars.output_precision
    );
    
 //   time(&solver.time_finish);
    if (verbose == ON) printf("Finished.\n\n");
    
}

