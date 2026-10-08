/*
*****************************************************************************
LOAD PARAMETER FILE PARAMETERS
---------------------

Load the filenames and options from the parameter file (.par) and defines the
simulation states. Detailed explanation of model parameters, files and simulation
states found in the manual. TJF.

*****************************************************************************
*/

#include "lisflood.h"

#define MODE_STRING(mode) (mode == CMD_LINE?"(Command Line)":"")

int ReadVerboseMode(int argc, char *argv[])
{
	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-v")) return ON;
	}

	return OFF;
}

void ReadConfiguration
(
	int argc,
	char *argv[],
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	const int verbose
)
{
	// Assume the last command line argument is the parameter file, check it doesn't start with a '-'
	if (argc == 1 || argv[argc - 1][0] == '-')
	{
		printf("no parameter file specified\n");
		exit(-1);
	}

	ReadParamFile(argv[argc-1], Fnameptr, Statesptr, Parptr, Solverptr, verbose);
	if (verbose == ON) printf("Reading parameters done.\n\n");
	if (verbose == ON) printf("Reading command line.\n");
	ReadCommandLine(argc, argv, Fnameptr, Statesptr, Parptr, Solverptr, verbose);
	if (verbose == ON) printf("Reading command line done.\n\n");

	CheckParams(Fnameptr, Statesptr, Parptr, Solverptr, verbose);
}

int read_string_param(const char * current_param_name, const int line_number, char* param_value_ptr, const char* param_name, char* filename_ptr, int verbose, int mode)
{
	if (STRCMPi(param_name, current_param_name) == 0)
	{
		/*int ret = */sscanf(param_value_ptr, "%255s", filename_ptr);
		//if (ret != 1)
		//{
		//	printf("error reading filename param %s param_value_ptr: %d\n", current_param_name, line_number);
		//	exit(-1);
		//}
		if (strlen(filename_ptr) == 255)
		{
			printf("file name too long %s line: %d max (255)\n", current_param_name, line_number);
			exit(-1);
		}
		if (verbose == ON) printf("%s = %s\n", param_name, filename_ptr);
		return 1;
	}
	return 0;
}

int read_numeric_param(const char * current_param_name, const int line_number, char* param_value_ptr, const char* param_name, NUMERIC_TYPE* numeric_param_ptr, int verbose, int mode)
{
	if (STRCMPi(param_name, current_param_name) == 0)
	{
		int ret = sscanf(param_value_ptr, "%" NUM_FMT"", numeric_param_ptr);
		if (ret != 1)
		{
			printf("error reading decimal param %s line: %d\n", current_param_name, line_number);
			exit(-1);
		}
		if (verbose == ON) printf("%s = %.3" NUM_FMT"\n", param_name, *numeric_param_ptr);
		return 1;
	}
	return 0;
}

/// reads option [value] or option
/// tries to read the option value, otherwise used default_value
int read_numeric_param(const char * current_param_name, const int line_number, char* param_value_ptr, const char* param_name, NUMERIC_TYPE* numeric_param_ptr, int verbose, int mode, NUMERIC_TYPE default_value)
{
	if (STRCMPi(param_name, current_param_name) == 0)
	{
		int ret = sscanf(param_value_ptr, "%" NUM_FMT"", numeric_param_ptr);
		if (ret != 1)
		{
			(*numeric_param_ptr) = default_value;
		}
		if (verbose == ON) printf("%s = %.3" NUM_FMT"\n", param_name, *numeric_param_ptr);
		return 1;
	}
	return 0;
}

int read_integer_param(const char * current_param_name, const int line_number, char* param_value_ptr, const char* param_name, int* integer_param_ptr, int verbose, int mode)
{
	if (STRCMPi(param_name, current_param_name) == 0)
	{
		int ret = sscanf(param_value_ptr, "%i", integer_param_ptr);
		if (ret != 1)
		{
			printf("error reading integer param %s line: %d\n", current_param_name, line_number);
			exit(-1);
		}
		if (verbose == ON) printf("%s = %d\n", param_name, *integer_param_ptr);
		return 1;
	}
	return 0;
}

int read_empty_param(const char * current_param_name, const int line_number, char* param_value_ptr, const char* param_name, int verbose, int mode)
{
	if (STRCMPi(param_name, current_param_name) == 0)
	{
		if (verbose == ON) printf("%s\n", param_name);
		return 1;
	}
	return 0;
}

//mode PARAM_FILE read from file
//mode PARAM_CMD read from command line
void CheckParam(char* param_name, char* param_value_ptr, int line_number, Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose, int mode)
{
	//ignore these params - already handled
	if (STRCMPi("v", param_name) == 0 ||
		STRCMPi("nch", param_name) == 0)
		return;

	if (read_string_param(param_name, line_number, param_value_ptr, "DEMfile", Fnameptr->demfilename, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "startfile", Fnameptr->startfilename, verbose, mode))
	{
		Statesptr->startfile = ON;
		Statesptr->binarystartfile = OFF;
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "binarystartfile", Fnameptr->startfilename, verbose, mode))
	{
		Statesptr->startfile = OFF;
		Statesptr->binarystartfile = ON;
		if (verbose == ON) printf("Using binary startfile\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "startelev", verbose, mode))
	{
		Statesptr->startelev = ON;
		if (verbose == ON) printf("Using water surface elevations for restart\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "resroot", Fnameptr->res_prefix, verbose, mode))
	{
		if (verbose == ON && mode == CMD_LINE) printf("Results root name reset by command line: %s\n", Fnameptr->res_prefix);
		return;
	}

	if (read_string_param(param_name, line_number, param_value_ptr, "dirroot", Fnameptr->res_dirname, verbose, mode) ||
		read_string_param(param_name, line_number, param_value_ptr, "dir", Fnameptr->res_dirname, verbose, mode))
	{
		if (verbose == ON && mode == CMD_LINE) printf("Output folder set by command line: %s\n", Fnameptr->res_dirname);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "fpfric", &Parptr->FPn, verbose, mode) ||
		read_numeric_param(param_name, line_number, param_value_ptr, "nfp", &Parptr->FPn, verbose, mode))
	{
		if (verbose == ON && mode == CMD_LINE) printf("Floodplain friction reset by command line: %" NUM_FMT"\n", Parptr->FPn);
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "startq", verbose, mode))
	{
		Statesptr->startq = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "startq2d", verbose, mode))
	{
		Statesptr->startq2d = ON;
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "ch_start_h", &Parptr->ch_start_h, verbose, mode))
	{
		Statesptr->start_ch_h = ON;
		return;
	}

	if (read_numeric_param(param_name, line_number, param_value_ptr, "sim_time", &Solverptr->Sim_Time, verbose, mode) ||
		read_numeric_param(param_name, line_number, param_value_ptr, "simtime", &Solverptr->Sim_Time, verbose, mode))
		return;

	if (read_numeric_param(param_name, line_number, param_value_ptr, "initial_tstep", &Solverptr->InitTstep, verbose, mode))
		return;
	if (read_numeric_param(param_name, line_number, param_value_ptr, "massint", &Parptr->MassInt, verbose, mode))
		return;
	if (read_numeric_param(param_name, line_number, param_value_ptr, "saveint", &Parptr->SaveInt, verbose, mode))
		return;
	if (read_numeric_param(param_name, line_number, param_value_ptr, "htol", &Solverptr->htol, verbose, mode))
		return;
	if (read_numeric_param(param_name, line_number, param_value_ptr, "qlimfact", &Solverptr->Qlimfact, verbose, mode))
	{
		if (verbose == ON) printf("Relaxing Q limit by factor x%" NUM_FMT" \n\n", Solverptr->Qlimfact);
		return;
	}

	if (read_integer_param(param_name, line_number, param_value_ptr, "ts_multiple", &Solverptr->ts_multiple, verbose, mode))
	{
		if (verbose == ON) printf("Running channel at x%d multiples of floodplain timestep\n\n", Solverptr->ts_multiple);
		return;
	}

	if (read_numeric_param(param_name, line_number, param_value_ptr, "overpass", &Parptr->op, verbose, mode))
	{
		Statesptr->single_op = ON;
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "overpassfile", Fnameptr->opfilename, verbose, mode))
	{
		Statesptr->multi_op = ON;
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "Qfile", Fnameptr->qfilename, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "manningfile", Fnameptr->nfilename, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCmanningfile", Fnameptr->SGCnfilename, verbose, mode))
	{
		if (verbose == ON) printf("SGC Manning file %s\n", Fnameptr->SGCnfilename);
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCdirnfile", Fnameptr->SGCdirnfilename, verbose, mode))
	{
		if (verbose == ON) printf("SGC flow directions file %s\n", Fnameptr->SGCnfilename);
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "riverfile", Fnameptr->rivername, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "bcifile", Fnameptr->bcifilename, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "bdyfile", Fnameptr->bdyfilename, verbose, mode))
		return;
	if (read_string_param(param_name, line_number, param_value_ptr, "weirfile", Fnameptr->weirfilename, verbose, mode) ||
		read_string_param(param_name, line_number, param_value_ptr, "weir", Fnameptr->weirfilename, verbose, mode))
	{
		if (verbose == ON && mode == CMD_LINE) printf("Weir file set by command line: %s\n", Fnameptr->weirfilename);
		return;
	}
	// Turn on steady-state checking, and check for a specified tolerance (MDW)
	if (read_empty_param(param_name, line_number, param_value_ptr, "steady", verbose, mode))
	{
		Statesptr->steadycheck = ON;
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "steadytol", &Parptr->steadyQtol, verbose, mode))
	{
		return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "ch_dynamic", verbose, mode) ||
		read_empty_param(param_name, line_number, param_value_ptr, "dynsw", verbose, mode))
	{
		Solverptr->dynsw = 1;
		if (verbose == ON) printf("Full dynamic steady state channel solver turned on for initial solution\n");
		return;
	}

	if (read_string_param(param_name, line_number, param_value_ptr, "porfile", Fnameptr->porfilename, verbose, mode))
	{
		Statesptr->porosity = ON;
		return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "depthoff", verbose, mode))
	{
		Statesptr->save_depth = OFF;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "elevoff", verbose, mode))
	{
		Statesptr->save_elev = OFF;
		return;
	}
	
	if (read_empty_param(param_name, line_number, param_value_ptr, "vtkoff", verbose, mode))
	{
		Statesptr->save_vtk = OFF;
		return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "adaptoff", verbose, mode))
	{
		Statesptr->adaptive_ts = OFF;
		Statesptr->qlim = ON;
		if (verbose == ON) printf("Using Qlim formulation for floodplain flow\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "acceleration", verbose, mode))
	{
		Statesptr->acceleration = ON;
		Statesptr->adaptive_ts = OFF;
		Statesptr->qlim = OFF;
		if (verbose == ON)	printf("\nUsing acceleration formulation for floodplain flow\n");

		return;
	}

	// Don't combine adaptoff and acceleration keywords. Program aborts.
	if ((read_empty_param(param_name, line_number, param_value_ptr, "acceleration", verbose, mode) && Statesptr->qlim == ON) ||
		(read_empty_param(param_name, line_number, param_value_ptr, "adaptoff", verbose, mode) && Statesptr->acceleration == ON))
	{
		if (verbose == ON) printf("\nIncompatible timestep functions chosen (adaptoff and acceleration) - Check parameter file. Aborting\n");
		exit(1);
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "chainageoff", verbose, mode))
	{
		Statesptr->chainagecalc = OFF;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "qoutput", verbose, mode)) {
		Statesptr->save_Qs = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "maxdepthonly", verbose, mode))
	{
		Statesptr->maxdepthonly = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "voutput", verbose, mode))
	{
		Statesptr->voutput = ON;
		Statesptr->voutput_stage = ON; // only applicable if stages already configured
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "voutput_max", verbose, mode))
	{
		Statesptr->voutput = ON;
		Statesptr->voutput_max = ON; // max only written if explicitly enabled (saves having to calculate velocity at each time step
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "voutput_stage", verbose, mode))
	{
		Statesptr->voutput_stage = ON;// only applicable if stages already configured - option to save velocity at stage without having to enable voutput
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "hazard", verbose, mode))
	{
		Statesptr->voutput = ON;
		Statesptr->voutput_max = ON;
		Statesptr->voutput_stage = ON; // only applicable if stages already configured
		Statesptr->hazard = ON;
		if (verbose == ON) printf("\nHazard mode on\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "profiles", verbose, mode))
	{
		Statesptr->profileoutput = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "debug", verbose, mode))
	{
		Statesptr->debugmode = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "comp_out", verbose, mode))
	{
		Statesptr->comp_out = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "standard_extensions_out", verbose, mode))
	{
		Statesptr->output_params.standard_extensions = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "netcdf_out", verbose, mode))
	{
		Statesptr->output_params.netcdf_out = ON;
		printf("Using binary raster output\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "binary_out", verbose, mode))
	{
		Statesptr->binary_out = ON;
		Statesptr->output_params.binary_out = ON;
		printf("Using binary raster output\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "ascii_out", verbose, mode))
	{
		Statesptr->output_params.ascii_out = ON;
		printf("Using ascii raster output\n");
		return;
	}
	// JN: added to request maxH, maxHtm totalHtm and initHtm be calulated at the mass interval
	if (read_empty_param(param_name, line_number, param_value_ptr, "mint_hk", verbose, mode))
	{
		Statesptr->mint_hk = ON;
		return;
	}
	// JN/IV: added to turn on Roe solever, don't combine with adaptoff
	if (read_empty_param(param_name, line_number, param_value_ptr, "Roe", verbose, mode))
	{
		Statesptr->Roe = ON;
		Statesptr->Roe_slow = OFF;
		Statesptr->adaptive_ts = OFF;
		Statesptr->acceleration = OFF;
		Statesptr->qlim = OFF;
		if (verbose == ON) printf("\nUsing Roe formulation for floodplain flow\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "Roe_slow", verbose, mode))
	{
		Statesptr->Roe_slow = ON; // uses ghost cell version of Roe solver.
		Statesptr->Roe = ON;
		Statesptr->adaptive_ts = OFF;
		Statesptr->acceleration = OFF;
		Statesptr->qlim = OFF;
		if (verbose == ON) printf("\nUsing Roe formulation for floodplain flow\n");
		return;
	}
	if ((read_empty_param(param_name, line_number, param_value_ptr, "Roe", verbose, mode) && Statesptr->qlim == ON) ||
		(read_empty_param(param_name, line_number, param_value_ptr, "adaptoff", verbose, mode) && Statesptr->Roe == ON) ||
		(read_empty_param(param_name, line_number, param_value_ptr, "acceleration", verbose, mode) && Statesptr->Roe == ON))
	{
		if (verbose == ON) printf("\nIncompatible timestep functions chosen (adaptoff/qlim/acceleration and Roe) - Check parameter file. Aborting\n");
		exit(1);
	}
	// MT: added to output adaptive timestep & qlimits & diffusive channel solver
	if (read_empty_param(param_name, line_number, param_value_ptr, "toutput", verbose, mode))
	{
		Statesptr->save_Ts = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "qloutput", verbose, mode))
	{
		Statesptr->save_QLs = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "diffusive", verbose, mode))
	{
		Statesptr->diffusive = ON;
		if (verbose == ON) printf("\nUsing diffusive solver for channel flow\n");
		return;
	}

	if (read_string_param(param_name, line_number, param_value_ptr, "ascheader", Fnameptr->ascheaderfilename, verbose, mode))
	{
		Statesptr->alt_ascheader = ON;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "calcarea", verbose, mode))
	{
		Statesptr->calc_area = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "calcmeandepth", verbose, mode))
	{
		Statesptr->calc_meandepth = ON;
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "calcvolume", verbose, mode))
	{
		Statesptr->calc_volume = ON;
		return;
	}

	if (read_string_param(param_name, line_number, param_value_ptr, "stagefile", Fnameptr->stagefilename, verbose, mode))
	{
		Statesptr->save_stages = ON;
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "gaugefile", Fnameptr->gaugefilename, verbose, mode))
	{
		Statesptr->gsection = ON;
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "infiltration", &Parptr->InfilRate, verbose, mode) ||
		read_numeric_param(param_name, line_number, param_value_ptr, "inf", &Parptr->InfilRate, verbose, mode))
	{
		Statesptr->calc_infiltration = ON;
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "infilfile", Fnameptr->infilfilename, verbose, mode))
	{
		Statesptr->calc_infiltration = ON;
		Statesptr->calc_distributed_infiltration = ON;
		return;
	}

	if (read_string_param(param_name, line_number, param_value_ptr, "evaporation", Fnameptr->evapfilename, verbose, mode))
	{
		Statesptr->calc_evap = ON;
		return;
	}

	// Enable rainfall CCS
	if (read_string_param(param_name, line_number, param_value_ptr, "rainfall", Fnameptr->rainfilename, verbose, mode))
	{
		Statesptr->rainfall = ON;
		return;
	}

	//Enable distributed rainfall mask
	if (read_string_param(param_name, line_number, param_value_ptr, "rainfallmask", Fnameptr->rainmaskname, verbose, mode))
	{
		Statesptr->rainfallmask = ON;
		return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "routing", verbose, mode)) // Enable routing scheme CCS
	{
		Statesptr->routing = ON;
		return;
	}
	// Enable routing scheme mass balance check //Toby Dunne
	// Checks each routing cell, and if total volume out exceeds total volume, reduces the routing Q
	// Enabling this check slows the model processing ~20%
	if (read_empty_param(param_name, line_number, param_value_ptr, "routing_mass_check", verbose, mode))
	{
		Statesptr->routing = ON;
		Statesptr->routing_mass_check = ON;
		return;
	}

	// Legacy routing keyword CCS
	if (read_string_param(param_name, line_number, param_value_ptr, "rainfallrouting", Fnameptr->rainfilename, verbose, mode))
	{
		Statesptr->rainfall = ON;
		Statesptr->routing = ON;
		return;
	}

	// Set speed at which shallow water (<depththresh) is routed across DEM CCS 14/03/2012
	if (read_numeric_param(param_name, line_number, param_value_ptr, "routingspeed", &Parptr->Routing_Speed, verbose, mode))
	{
		Statesptr->routing = ON;
		return;
	}

	// Set slope where routing replaces shallow water eqn  CCS 14/03/2012
	if (read_numeric_param(param_name, line_number, param_value_ptr, "routesfthresh", &Parptr->RouteSfThresh, verbose, mode))
	{
		Statesptr->routing = ON;
		return;
	}

	// Set slope where routing replaces shallow water eqn  CCS 14/03/2012
	if (read_numeric_param(param_name, line_number, param_value_ptr, "diffusive_froude_thresh", &Parptr->DiffusiveFroudeThresh, verbose, mode))
	{
		Statesptr->diffusive_switch = ON;
		return;
	}

	if (read_numeric_param(param_name, line_number, param_value_ptr, "checkpoint", &Parptr->checkfreq, verbose, mode, C(1.0)))
	{
		Statesptr->checkpoint = ON;
		if (Parptr->checkfreq <= C(0.0)) Parptr->checkfreq = CHKINTERVAL; // if interval is less than zero set to default
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "loadcheck", Fnameptr->loadCheckpointFilename, verbose, mode))
	{
		Statesptr->checkpoint = ON;
		Statesptr->checkfile = ON;
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "resettimeinit", &Parptr->reset_timeinit_time, verbose, mode))
	{
		Statesptr->reset_timeinit = ON;
		if (verbose == ON) printf("\n Time of initial inundation will be reset at %" NUM_FMT" seconds\n", Parptr->reset_timeinit_time);
		return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "fv1", verbose, mode))
    {
		Statesptr->fv1 = ON;
        Solverptr->cfl = C(0.5);
        Solverptr->DepthThresh = C(1e-3);
		if (verbose == ON) printf("Using FV1 2D SWE solver\n");
        return;
    }
	else if (read_empty_param(param_name, line_number, param_value_ptr, "fv2", verbose, mode))
	{
		Statesptr->fv2 = ON;
		Solverptr->cfl = C(0.33);
		Solverptr->DepthThresh = C(1e-3);
		//Solverptr->FV2DepthThresh = C(-1.0);
		//Solverptr->FV2ThinDepthTstep = C(0.0);
		if (verbose == ON) printf("Using FV2 2D SWE solver\n");
		return;
	}	
	else if (read_empty_param(param_name, line_number, param_value_ptr, "acc_nugrid", verbose, mode))
	{
		Statesptr->acc_nugrid = ON;
		Solverptr->cfl = C(0.7);
		Solverptr->DepthThresh = C(1e-3);
		//Solverptr->FV2DepthThresh = C(-1.0);
		//Solverptr->FV2ThinDepthTstep = C(0.0);
		if (verbose == ON) printf("\nUsing acceleration formulation on non-uniform grid\n");
		return;
	}
	else if (read_empty_param(param_name, line_number, param_value_ptr, "mwdg2", verbose, mode))
	{
		Statesptr->mwdg2 = ON;
		Solverptr->cfl = C(0.3);
		Solverptr->DepthThresh = C(1e-3);
		//Solverptr->FV2DepthThresh = C(-1.0);
		//Solverptr->FV2ThinDepthTstep = C(0.0);
		if (verbose == ON) printf("\nUsing adaptive MWDG2 2D SWE solver\n");
		return;
	}
	else if (read_empty_param(param_name, line_number, param_value_ptr, "hwfv1", verbose, mode))
	{
		Statesptr->hwfv1 = ON;
		Solverptr->cfl = C(0.5);
		Solverptr->DepthThresh = C(1e-3);
		//Solverptr->FV2DepthThresh = C(-1.0);
		//Solverptr->FV2ThinDepthTstep = C(0.0);
		if (verbose == ON) printf("\nUsing adaptive HWFV1 2D SWE solver\n");
		return;
	}
	else if (read_empty_param(param_name, line_number, param_value_ptr, "dg2", verbose, mode))
	{
		Statesptr->dg2 = ON;
        Solverptr->cfl = C(0.33);
        Solverptr->DepthThresh = C(1e-3);
        Solverptr->DG2DepthThresh = C(-1.0);
        Solverptr->DG2ThinDepthTstep = C(0.0);
		if (verbose == ON) printf("Using DG2 2D SWE solver\n");
        return;
	}

	// Reset the cfl condition for acceleration version - reduce to increase stability
	if (read_numeric_param(param_name, line_number, param_value_ptr, "cfl", &Solverptr->cfl, verbose, mode))
	{
		if (mode == CMD_LINE && Statesptr->acceleration == OFF)
			printf("WARNING: CFL value changed on command line but acceleration version off\n");
		if (verbose == ON) printf("cfl changed to %" NUM_FMT" %s\n", Solverptr->cfl, MODE_STRING(mode));
		return;
	}

	// Reset the theta parameter for acceleration version - reduce to increase numerical diffusion
	if (read_numeric_param(param_name, line_number, param_value_ptr, "theta", &Solverptr->theta, verbose, mode))
	{
		if (verbose == ON)
		{
			if (mode == CMD_LINE && Statesptr->acceleration == OFF)
				printf("WARNING: theta value changed on command line but acceleration version off\n");
			printf("theta changed to %" NUM_FMT" %s\n", Solverptr->theta, MODE_STRING(mode));
		}
		return;
	}

	// Uses the 1D (old) version of the friction term
	if (read_empty_param(param_name, line_number, param_value_ptr, "1Dfriction", verbose, mode))
	{
		Solverptr->fricSolver2D = OFF;
		if (verbose == ON) printf("Using the 1D version of the friction term\n");
		return;
	}

	// Reset the dhlin condition for adaptive version - reduce to increase stability/ increase to reduce computation time
	// Overwrites the standard value which is set as gradient = C(0.0002) (i.e. dhlin is a function of dx)
	if (read_numeric_param(param_name, line_number, param_value_ptr, "dhlin", &Solverptr->dhlin, verbose, mode))
	{
		Statesptr->dhoverw = ON;
		if (verbose == ON && mode == CMD_LINE) printf("dhlin value reset by command line: %" NUM_FMT"\n", Solverptr->dhlin);
		return;
	}
	// Options to change depth and momentum thresholds defaults used if not set
	if (read_numeric_param(param_name, line_number, param_value_ptr, "depththresh", &Solverptr->DepthThresh, verbose, mode))
	{
		if (verbose == ON) printf("Depth threshold changed to %" NUM_FMT" \n", Solverptr->DepthThresh);
		return;
	}

	// adaptation
	// Options to change adaptation threshold defaults used if not set
	if (read_numeric_param(param_name, line_number, param_value_ptr, "epsilon", &Solverptr->epsilon, verbose, mode))
	{
		if (verbose == ON) printf("Adaptation threshold changed to %" NUM_FMT" \n", Solverptr->epsilon);
		return;
	}

	// Options to change maximum resolution level defaults used if not set
	if (read_integer_param(param_name, line_number, param_value_ptr, "L", &Solverptr->L, verbose, mode))
	{
		if (verbose == ON) printf("Maximum resolution level changed to %" NUM_FMT" \n", Solverptr->L);
		return;
	}









	if (read_numeric_param(param_name, line_number, param_value_ptr, "dg2depththresh", &Solverptr->DG2DepthThresh, verbose, mode))
	{
        if (Statesptr->dg2 == OFF)
        {
            printf("\ndg2depththresh only applicable to dg2 solver. Aborting...\n\n");
            exit(1);
        }
		if (verbose == ON) printf("DG2 depth threshold changed to %" NUM_FMT" \n", Solverptr->DG2DepthThresh);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "dg2thindepth_tstep", &Solverptr->DG2ThinDepthTstep, verbose, mode))
	{
        if (Statesptr->dg2 == OFF)
        {
            printf("\ndg2thindepth_tstep only applicable to dg2 solver. Aborting...\n\n");
            exit(1);
        }
		if (verbose == ON) printf("DG2 thin depth timestep changed to %" NUM_FMT" \n", Solverptr->DG2ThinDepthTstep);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "momentumthresh", &Solverptr->MomentumThresh, verbose, mode))
	{
		if (verbose == ON) printf("Momentum threshold changed to %" NUM_FMT" \n", Solverptr->MomentumThresh);
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "drycheckoff", verbose, mode))
	{ // turns DryCheck off use at own risk!
		Statesptr->drychecking = OFF;
		if (verbose == ON) printf("DryCheck is off\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "multiriverfile", Fnameptr->multiriverfilename, verbose, mode))
	{
		Statesptr->multiplerivers = ON;
		if (verbose == ON) printf("\nMultiple river mode selected\n");
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCp", &Parptr->SGC_p, verbose, mode))
	{
		Statesptr->SGC = ON;
		if (verbose == ON) printf("SGC exponent changed to %.5" NUM_FMT" \n", Parptr->SGC_p);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCr", &Parptr->SGC_r, verbose, mode))
	{
		Statesptr->SGC = ON;
		if (verbose == ON) printf("SGC multiplier changed to %.5" NUM_FMT" \n", Parptr->SGC_r);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCm", &Parptr->SGC_m, verbose, mode))
	{
		if (verbose == ON) printf("SGC meander coefficient changed to %.5" NUM_FMT" \n", Parptr->SGC_m);
		return;
	}
	if (read_integer_param(param_name, line_number, param_value_ptr, "SGCchan", &Parptr->SGCchan_type, verbose, mode))
	{
		if (verbose == ON) printf("SGC channel changed to %i \n", Parptr->SGCchan_type);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCs", &Parptr->SGC_s, verbose, mode))
	{
		if (verbose == ON) printf("SGCs changed to %.5" NUM_FMT" \n", Parptr->SGC_s);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGC2", &Parptr->SGC_2, verbose, mode))
	{
		if (verbose == ON) printf("SGC2 changed to %.5" NUM_FMT" \n", Parptr->SGC_2);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCn", &Parptr->SGC_n, verbose, mode))
	{
		if (verbose == ON) printf("SGCn changed to %.5" NUM_FMT" \n", Parptr->SGC_n);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "SGCa", &Parptr->SGC_a, verbose, mode))
	{
		if (verbose == ON) printf("SGCa changed to %.5" NUM_FMT" \n", Parptr->SGC_a);
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "SGCbfh_mode", verbose, mode))
	{
		Statesptr->SGCbfh_mode = ON; // Turns on mode where parameter p is used as the channel bank full depth
		if (verbose == ON) printf("Using sub-grid parameter p as bank full depth\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "SGCA_mode", verbose, mode))
	{
		Statesptr->SGCbfh_mode = ON; // Turns on mode where parameter p is used as the channel bank full area
		if (verbose == ON) printf("Using sub-grid parameter p as bank full area\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "SGCvoutput", verbose, mode))
	{
		Statesptr->SGCvoutput = ON; // Turns on mode where parameter p is used as the channel bank full area
		if (verbose == ON) printf("Output sub-grid channel velocity\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "SGC_enable", verbose, mode))
	{
		// enables the sgm_fast optimized version, without necessarily having any channels
		Statesptr->SGC = ON;
		Statesptr->acceleration = ON; // sug-grid channels only work with acceleration model so ensure its in this model.
		Statesptr->adaptive_ts = OFF;
		Statesptr->qlim = OFF;
		if (verbose == ON) printf("Using sub-grid channels and acceleration formulation\n");
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCwidth", Fnameptr->SGCwidthfilename, verbose, mode))
	{
		Statesptr->SGC = ON;
		Statesptr->acceleration = ON; // sug-grid channels only work with acceleration model so ensure its in this model.
		Statesptr->adaptive_ts = OFF;
		Statesptr->qlim = OFF;
		if (verbose == ON) printf("Using sub-grid channels and acceleration formulation\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCbed", Fnameptr->SGCbedfilename, verbose, mode))
	{
		Statesptr->SGCbed = ON;
		if (verbose == ON) printf("SGC bed elevation read\n");
		return;
	}
	//if (read_string_param(param_name, line_number, param_value_ptr, "SGClevee", Fnameptr->SGCleveefilename, verbose, mode))
	//{
	//	Statesptr->SGClevee = ON;
	//	if (verbose == ON) printf("SGC levee read\n");
	//	return;
	//}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCcat_area", Fnameptr->SGCcat_areafilename, verbose, mode))
	{
		Statesptr->SGCcat_area = ON;
		if (verbose == ON) printf("SGC catchment area read\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCchangroup", Fnameptr->SGCchangroupfilename, verbose, mode))
	{
		Statesptr->SGCchangroup = ON;
		if (verbose == ON) printf("SGC channel type\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCchanprams", Fnameptr->SGCchanpramsfilename, verbose, mode))
	{
		Statesptr->SGCchanprams = ON;
		if (verbose == ON) printf("SGC channel parameters\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "SGCbank", Fnameptr->SGCbankfilename, verbose, mode))
		return;
	if (read_numeric_param(param_name, line_number, param_value_ptr, "tstart", &Solverptr->t, verbose, mode))
	{
		Parptr->SaveTotal = Solverptr->t;
		Parptr->MassTotal = Solverptr->t;
		Parptr->maxintTotal = Solverptr->t;
		if (verbose == ON) printf("Simulation start time changed to %" NUM_FMT" \n", Solverptr->t);
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "gravity", &Solverptr->g, verbose, mode))
	{
		if (verbose == ON) printf("Simulation gravity changed to %" NUM_FMT" \n", Solverptr->g);
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "latlong", verbose, mode))
	{
		Statesptr->latlong = ON;
		if (verbose == ON) printf("Lat-Long coordinate system on.\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "SGCd8", verbose, mode))
	{
		Statesptr->SGCd8 = ON;
		if (verbose == ON) printf("Sub Grid Channel D8 directions on.\n");
		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "dist_routing", verbose, mode))
	{
		Statesptr->dist_routing = ON;
		if (verbose == ON) printf("Using slope dependent routing velocity.\n");
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "kill", &Parptr->killsim_time, verbose, mode))
	{
		if (verbose == ON) printf("\n***** WARNING: Simulation will be killed after %.2" NUM_FMT" hours *****\n", Parptr->killsim_time);
		Parptr->killsim_time *= 3600;//convert hours to seconds
		Statesptr->killsim = ON;
		return;
	}

	//find out if we are going to compress output on the fly
	if (read_empty_param(param_name, line_number, param_value_ptr, "gzip", verbose, mode))
	{
		Statesptr->call_gzip = ON;
		Statesptr->output_params.call_gzip = ON;
		if (verbose == ON) printf("\nOutput will be compressed using Gzip\n");
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "log", verbose, mode))
	{
		Statesptr->logfile = ON;
		return;
	}
	
	if (read_string_param(param_name, line_number, param_value_ptr, "damfile", Fnameptr->Damparfilename, verbose, mode))
	{
		if (verbose == ON)
			if (mode == CMD_LINE) printf("Dam file set by command line: %s\n", Fnameptr->Damparfilename);
			else printf("Dam file set by par file: %s\n", Fnameptr->Damparfilename);
			Statesptr->DamMode = ON;
			return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "dammask", Fnameptr->DamMaskfilename, verbose, mode))
	{
		if (verbose == ON)
			if (mode == CMD_LINE) printf("Dam mask set by command line: %s\n", Fnameptr->DamMaskfilename);
			else printf("Dam mask set by par file: %s\n", Fnameptr->DamMaskfilename);
			Statesptr->DammaskRead = ON;
			return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "max_Froude", &Parptr->max_Froude, verbose, mode))
	{
		if (verbose == ON) printf("Local inertia max V changed to Froude %.5" NUM_FMT" \n", Parptr->max_Froude);
#if _CALCULATE_Q_MODE != 1
		printf("ERROR Local inertia max V not enbled in this exectable CALCULATE_Q_MODE=1 needed during compilation\n");
		printf("To run this model with current executable remove max_Froude from parfile\n");
		exit(EXIT_FAILURE);
#endif

		return;
	}
	if (read_empty_param(param_name, line_number, param_value_ptr, "saveint_max", verbose, mode))
	{
		Statesptr->saveint_max = ON;
		Statesptr->maxint = OFF; 
		if (verbose == ON) printf("Output and reset max at save interval\n");
		return;
	}
	if (read_numeric_param(param_name, line_number, param_value_ptr, "maxint", &Parptr->maxint, verbose, mode))
	{
		Statesptr->saveint_max = OFF; 
		Statesptr->maxint = ON; 
		if (verbose == ON) printf("Output and reset max at maxint interval\n");
		return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "chanmask", Fnameptr->ChanMaskfilename, verbose, mode))
	{
		if (verbose == ON)
			if (mode == CMD_LINE) printf("Channel mask set by command line: %s\n", Fnameptr->ChanMaskfilename);
			else printf("Channel mask set by par file: %s\n", Fnameptr->ChanMaskfilename);
			Statesptr->ChanMaskRead = ON;
			return;
	}
	if (read_string_param(param_name, line_number, param_value_ptr, "linklist", Fnameptr->LinkListfilename, verbose, mode))
	{
		if (verbose == ON)
			if (mode == CMD_LINE) printf("Link list set by command line: %s\n", Fnameptr->LinkListfilename);
			else printf("link list set by par file: %s\n", Fnameptr->LinkListfilename);
			Statesptr->LinkListRead = ON;
			return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "cuda", verbose, mode))
	{
		Statesptr->cuda = ON;
		if (verbose == ON) printf("CUDA enabled\n");
        return;
	}

	if (read_empty_param(param_name, line_number, param_value_ptr, "limitslopes", verbose, mode))
    {
        if (Statesptr->dg2 == OFF && Statesptr->mwdg2 == OFF && Statesptr->fv2 == OFF)
        {
			printf("ERROR Slope limiter only applicable to DG2, MWDG2 and FV2 solvers\n");
			exit(EXIT_FAILURE);
        }
		Parptr->limit_slopes = ON;
		if (verbose == ON) printf("DG2/FV2 slope limiter enabled\n");
        return;
    }
	if (read_numeric_param(param_name, line_number, param_value_ptr, "krivodonovathresh", &Solverptr->krivodonova_threshold, verbose, mode))
	{
		if (verbose == ON) printf("Krivodonova threshold changed to %" NUM_FMT" \n", Solverptr->krivodonova_threshold);
		return;
	}

	if (read_integer_param(param_name, line_number, param_value_ptr, "output_precision", &Parptr->output_precision, verbose, mode))
	{
		if (verbose == ON) printf("Output precision changed to %i\n", Parptr->output_precision);
		return;
	}
    
	if (read_string_param(param_name, line_number, param_value_ptr, "dynamicrainfile", Fnameptr->dynamicrainfilename, verbose, mode))
	{
		Statesptr->dynamicrainfall = ON;
#if (_NETCDF != 1)			
		fprintf(stderr, "ERROR: dynamic rain requires lisflood to be compiled with NetCDF support\n");
        exit(1);
#endif
		return;
	}

    if (read_numeric_param(param_name, line_number, param_value_ptr, "nodata_elevation", &Parptr->nodata_elevation, verbose, mode))
    {
        if (verbose == ON) printf("DEM elevation set to %" NUM_FMT " for nodata values\n", Parptr->nodata_elevation);
        return;
    }

	if (read_empty_param(param_name, line_number, param_value_ptr, "drain_nodata", verbose, mode))
	{
		Parptr->drain_nodata = ON;
		if (verbose == ON) printf("Water will be removed from DEM NODATA cells\n");
        return;
	}

	// needs to be last in loop
	if (verbose == ON) printf("Unknown parameter ignored: %s.\n", param_name);
}

void ReadCommandLine(int argc, char *argv[], Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose)
{
	char* param_ptr;
	char* param_value;
	for (int i = 0; i < argc-1; i++)
	{
		if (argv[i][0] == '-')
		{
			param_ptr = &argv[i][1];
			if (i + 1 < argc)
				param_value = argv[i + 1];
			else
				param_value = "";
			CheckParam(param_ptr, param_value, i, Fnameptr, Statesptr, Parptr, Solverptr, verbose, CMD_LINE);
		}
	}
}

void CheckParams(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose)
{
	// do some basic checks to warn user of option conflicts
	if (Statesptr->output_params.ascii_out == OFF &&
		Statesptr->output_params.binary_out == OFF &&
		Statesptr->output_params.netcdf_out == OFF)
	{
		Statesptr->output_params.ascii_out = ON;
	}
	if (Statesptr->start_ch_h == ON && Statesptr->startq == ON && verbose == ON) printf("\nWARNING: startq option overides ch_start_h values\n\n");
	if (Statesptr->startfile == ON && Statesptr->startq == ON && verbose == ON) printf("\nWARNING: startfile H values overide startq values\n\n");
	if (Statesptr->startfile == ON && Statesptr->start_ch_h == ON && verbose == ON) printf("\nWARNING: startfile H values overide ch_start_h values\n\n");
	if (Statesptr->routing == ON && Statesptr->acceleration == OFF && Statesptr->SGC == OFF && verbose == ON) // CCS disable routing if not being used with inertial or SG model:
	{
		Statesptr->routing = OFF;
		printf("\nWARNING: Routing must be used with inertial or subgrid models. Routing disabled.\n\n");
	}
	if (Statesptr->latlong == ON && Statesptr->SGC == OFF && verbose == ON) // CCS abort if trying to use latlong without SG:
	{
		printf("\nWARNING: Latlong must be used with subgrid model. Aborting...\n\n");
		exit(1);
	}
	if (Statesptr->DamMode == ON && verbose == ON && Statesptr->SGC == OFF) // FEOL
	{
		printf("\nWARNING: Dam Mode Must be used with subgrid model. Aborting...\n\n");
		exit(1);
	}
	if (Statesptr->DamMode == ON && verbose == ON && Statesptr->DammaskRead == OFF) //FEOL
	{
		printf("\nWARNING: Dam Mask missing. Aborting...\n\n");
		exit(1);
	}

	if (Statesptr->DammaskRead == ON && Statesptr->DamMode == OFF && verbose == ON) //FEOL
	{
		printf("\nWARNING: Dam Params missing. Aborting ....\n\n");
		exit(1);
	}

#ifndef CUDA
	if (Statesptr->cuda == ON)
	{
		printf("\nERROR: lisflood has not been compiled with CUDA support. Aborting...\n\n");
		exit(1);
	}
#endif
}

//-----------------------------------------------------------------------------
// LOAD COMMAND LINE PARAMETERS INTO GLOBAL VARIABLES
void ReadParamFile(char *fname, Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose)
{
	char line_buffer[LINE_BUFFER_LEN];
	char* param_value_ptr;
	char param_name[80];
	FILE *par_fp;
	par_fp = fopen(fname, "r");

	if (par_fp == NULL)
	{
		printf("Parameter file '%s' not found. Aborting.\n", fname);
		fprintf(stderr, "Parameter file '%s' not found. Aborting.\n", fname);
		exit(-1);
	}

	if (verbose == ON) printf("Loading parameters... '%s'\n", fname);

	int line_number = 0;
	int ascii_enabled_by_param = OFF;
	while (fgets(line_buffer, LINE_BUFFER_LEN, par_fp))
	{
		line_number++;
		int line_len = strlen(line_buffer);
		if (line_len == LINE_BUFFER_LEN - 1)
		{
			printf("Parameter line too long line: %d\n", line_number);
			exit(-1);
		}
		if (line_len == 0 || line_buffer[0] == '#' || line_buffer[0] == '\n' ||
			(line_len > 1 && line_buffer[0] == '\r' && line_buffer[1] == '\n'))
			continue;

		//read the parameter name
		int ret = sscanf(line_buffer, "%79s ", param_name);
		if (ret != 1)
		{
			printf("Parameter name read error line: %d\n", line_number);
			exit(-1);
		}
		int param_name_len = strlen(param_name);
		if (param_name_len == 79)
		{
			printf("Parameter name too long line: %d\n", line_number);
			exit(-1);
		}
		param_value_ptr = line_buffer + param_name_len;

		CheckParam(param_name, param_value_ptr, line_number, Fnameptr, Statesptr, Parptr, Solverptr, verbose, PARAM_FILE);
	}

	fclose(par_fp);

	

	return;
}
