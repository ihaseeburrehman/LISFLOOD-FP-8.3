/*
#####################################################################################
LISFLOOD-FP flood inundation model
#####################################################################################

� copyright Bristol University Hydrology Research Group 2008

webpage -	http://www.ggy.bris.ac.uk/research/hydrology/models/lisflood
contact -	Professor Paul Bates, email: paul.bates@Bristol.ac.uk,
Tel: +44-117-928-9108, Fax: +44-117-928-7878

*/ 


#include "lisflood.h"
#include "VersionHistory.h"
#include "lisflood2/lisflood_processing.h"
#include "utility.h"
#include "sgc.h"
#include "swe/fv1.h"
#include "swe/dg2.h"
#ifdef CUDA
	#include "cuda/acc/cuda_acc_simulate.cuh"													
	#include "cuda/fv1/cuda_fv1_simulate.cuh"
	#include "cuda/dg2/cuda_dg2_simulate.cuh"
	#include "cuda/fv2/cuda_fv2_simulate.cuh"
    #include "cuda/acc_nugrid/cuda_acc_nugrid_simulate.cuh"
	#include "cuda/adaptive/cuda_adaptive_simulate.cuh"
#endif

#include "lisflood2/file_tool.h"

//---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
#ifdef TESTING
	RunTests();
	exit(0);
#endif

	int i, chseg;
	FILE *tmp_fp;
	char t1[255];
	NUMERIC_TYPE tmp;
	char tmp_sys_com[255]; // temporary string to hold system command

	// Instances of Structures
	Arrays Raster;
	Files Fps;
	Fnames ParFp;
	States SimStates;
	Pars Params;
	Solver ParSolver;
	BoundCs Bounds;
	Stage OutLocs;
	SGCprams SGCchanprams;
	DamData DamDataprams;


	//Instances of Vectors
	vector<ChannelSegmentType> ChannelSegments; // CCS Contains the channel information for ALL rivers (indexed using RiversIndex vector below).
	vector<ChannelSegmentType> *ChannelSegmentsVecPtr; // CCS
	ChannelSegmentsVecPtr = &ChannelSegments; // CCS

	vector<QID7_Store> QID7; //#CCS A temporary store for some ChannelSegments variables that cannot be written to the correct location during LoadRiver().
	vector<QID7_Store> *QID7_Vec_Ptr; // CCS
	QID7_Vec_Ptr = &QID7; // CCS

	vector<int> RiversIndex; // CCS Contains index values for ChannelSegments so we know where one river finishes and the next starts.
	vector<int> *RiversIndexVecPtr; // CCS
	RiversIndexVecPtr = &RiversIndex; // CCS

	memset(&Raster, 0, sizeof(Arrays));
	memset(&Fps, 0, sizeof(Files));
	memset(&ParFp, 0, sizeof(Fnames));
	memset(&SimStates, 0, sizeof(States));
	memset(&Params, 0, sizeof(Pars));
	memset(&ParSolver, 0, sizeof(Solver));
	memset(&Bounds, 0, sizeof(BoundCs));
	memset(&OutLocs, 0, sizeof(Stage));
	memset(&SGCchanprams, 0, sizeof(SGCprams));
	memset(&DamDataprams, 0, sizeof(DamData));

	// DEFINE & DECLARE: Pointers to Structures
	Arrays *Arrptr = &Raster;
	Fnames *Fnameptr = &ParFp;
	States *Statesptr = &SimStates;
	Pars *Parptr = &Params;
	Solver *Solverptr = &ParSolver;
	BoundCs *BCptr = &Bounds;
	Stage *Stageptr = &OutLocs;
	SGCprams *SGCptr = &SGCchanprams;
	DamData *Damptr = &DamDataprams; //FEOL

	// Define initial value for common simulation states (eg. verbose)
	SimStates.output_params = OutputParams();

	// Define initial value for parameters
	Params.dx = C(10.0);
	Params.dy = C(10.0);
	Params.dA = C(100.0);
	Params.tlx = C(0.0);
	Params.tly = C(0.0);
	Params.blx = C(0.0);
	Params.bly = C(0.0);
	Params.FPn = C(0.06);
	Params.SaveInt = C(1000.0);
	Params.SaveTotal = C(0.0);
	Params.MassInt = C(100.0);
	Params.MassTotal = C(0.0);
	Params.SaveNo = 0;
	Params.op = C(100.0);
	Params.InfilLoss = C(0.0);
	Params.EvapLoss = C(0.0);
	Params.RainLoss = C(0.0);
	Params.InfilTotalLoss = C(0.0);
	Params.EvapTotalLoss = C(0.0);
	Params.RainTotalLoss = C(0.0);
	Params.checkfreq = CHKINTERVAL;  // set default checkpointing file output interval
	Params.nextcheck = C(0.0);
	Params.reset_timeinit_time = 0;
	Params.op_multinum = 0; // default to zero or can cause problems with checkpointing if multipass not used
	Params.ch_start_h = C(2.0); // default start water depth for channel
	Params.steadyQtol = C(0.0005); // tolerance for steady-state definition
	Params.steadyInt = C(1800.0); // interval at which to assess steady-state
	Params.steadyTotal = C(0.0);
	Params.SGC_p = C(0.78); // default for sub grid channel exponent
	Params.SGC_r = C(0.12); // default for sub grid channel multiplier (British rivers average Hey and Thorne (1986))
	Params.SGCchan_type = 1; // defines the type of channel used by the sub_grid model, default is rectangular channel.
	Params.SGC_s = C(2.0); // sub-grid channel parameter used for some of the channel types, parabolic channel default.
	Params.SGC_2 = C(0.0); // sub-grid channel parameter used for some of the channel types, meaningless default.
	Params.SGC_n = C(0.035); // sub-grid channel parameter used for some of the channel types, meaningless default.
	Params.Routing_Speed = C(0.1); // CCS default routing speed for shallow rainfall flows C(0.1) m/s
	Params.RouteInt = C(0.0); // CCS will be reasigned when FlowDirDEM function is called
	Params.RouteSfThresh = C(0.1); // CCS water surface slope at which routing takes over from shallow water equations when routing is enabled.
	Params.SGC_m = 1;  // JCN meander coefficient for sub-grid model.
	Params.SGC_a = -1; // JCN upstream area for sub-grid model.
	Params.min_dx = C(10.0); // CCS Holds min_dx value (needed for variable dimension lat-long grids)
	Params.min_dy = C(10.0); // CCS Holds min_dy value (needed for variable dimension lat-long grids)
	Params.min_dx_dy = C(10.0); // CCS Holds min of min_dx and min_dy values (needed for variable dimension lat-long grids)
	Params.max_Froude = C(10000.0); // maximum Froude allowed in model, set way higher than will ever occure by default (JCN)
	Params.maxint = C(99999999999.0); // maxinterval save time
	Params.maxintTotal = C(0.0); // maxinterval save time
	Params.maxintcount = 0; // number of max ints
    Params.limit_slopes = OFF;
	Params.output_precision = 3;
    Params.nodata_elevation = DEM_NO_DATA;
    Params.drain_nodata = OFF;

	// Define initial values for boundary conditions
	Bounds.Qin = C(0.0);
	Bounds.Qout = C(0.0);
	Bounds.VolInMT = C(0.0);
	Bounds.VolOutMT = C(0.0);

	// Define initial values for arrays
	Raster.Manningsn = NULL;
	Raster.SGCManningsn = NULL;

	// Define initial values for solver settings
	ParSolver.Sim_Time = C(3600.0);
	ParSolver.InitTstep = C(10.0);		// Maximum timestep
	ParSolver.Nit = 360;
	ParSolver.itCount = 0;
	ParSolver.t = C(0.0);
	ParSolver.g = C(9.8065500000000);
	ParSolver.divg = (1 / (2 * ParSolver.g));
	ParSolver.cfl = C(0.7);
	ParSolver.SolverAccuracy = C(1e-4);
	ParSolver.dynsw = 0; // Switch for full dynamic steady state (1) or diffusive steady state (0)
	ParSolver.DepthThresh = C(1e-3);
	ParSolver.MomentumThresh = C(1e-2);
	ParSolver.MaxHflow = C(10.0);
	ParSolver.Hds = C(0.0);
	ParSolver.Qerror = C(0.0);
	ParSolver.Verror = C(0.0);
	ParSolver.dhlin = C(0.01);
	ParSolver.htol = C(1.0);
	ParSolver.Qlimfact = C(1.0);
	ParSolver.itrn_time = C(0.0);
	ParSolver.itrn_time_now = C(0.0);
	ParSolver.ts_multiple = 1;  // default to x1 timestep decouple multiple
	ParSolver.SGCtmpTstep = 1; // JCN any number 
	ParSolver.theta = C(1.0); // GAMA (for q-centred numerical scheme), C(1.0)= semi-implicit version (Bates et al 2010);
	ParSolver.fricSolver2D = ON; //GAMA: uses the 2D friction scheme as default
	ParSolver.maxH = C(0.0);
	ParSolver.krivodonova_threshold = C(10.0);
	ParSolver.SpeedThresh = C(1e-6);
	ParSolver.epsilon = C(0.001); // adaptation
	ParSolver.L = 0; // adaptation

	// Define default values for SimStates instance of States
	SimStates.diffusive = OFF;	// CCS added default state
	SimStates.ChannelPresent = OFF;
	SimStates.TribsPresent = ON;
	SimStates.NCFS = ON;
	SimStates.save_depth = ON;
	SimStates.save_elev = ON;
	SimStates.save_vtk = ON;
	SimStates.single_op = OFF;
	SimStates.multi_op = OFF;
	SimStates.calc_area = OFF;
	SimStates.calc_meandepth = OFF;
	SimStates.calc_volume = OFF;
	SimStates.save_stages = OFF;
	SimStates.adaptive_ts = ON;
	SimStates.qlim = OFF; //TJF: Switch for qlim version, default is OFF
	SimStates.acceleration = OFF; //PB: Switch for acceleration version, default is OFF
	SimStates.debugmode = OFF;
	SimStates.save_Qs = OFF;
	SimStates.calc_infiltration = OFF;
	SimStates.calc_distributed_infiltration = OFF;
	SimStates.call_gzip = OFF;
	SimStates.alt_ascheader = OFF;
	SimStates.checkpoint = OFF;
	SimStates.checkfile = OFF;
	SimStates.calc_evap = OFF;
	SimStates.rainfall = OFF;
	SimStates.rainfallmask = OFF;
	SimStates.routing = OFF; //CCS: Switch for rainfall routing routine 
	SimStates.reset_timeinit = OFF;
	SimStates.profileoutput = OFF;
	SimStates.porosity = OFF;
	SimStates.weirs = OFF;
	SimStates.save_Ts = OFF;
	SimStates.save_QLs = OFF;
	SimStates.startq = OFF;
	SimStates.logfile = OFF;
	SimStates.startfile = OFF;
	SimStates.start_ch_h = OFF;
	SimStates.comp_out = OFF;
	SimStates.chainagecalc = ON;
	SimStates.mint_hk = OFF;
	SimStates.Roe = OFF;
	SimStates.killsim = OFF;
	SimStates.dhoverw = OFF;
	SimStates.drychecking = ON;
	SimStates.voutput = OFF;
	SimStates.steadycheck = OFF;
	SimStates.hazard = OFF;
	SimStates.startq2d = OFF;
	SimStates.Roe_slow = OFF;
	SimStates.multiplerivers = OFF;
	SimStates.SGC = OFF;
	SimStates.SGCbed = OFF;
	SimStates.SGClevee = OFF;
	SimStates.SGCcat_area = OFF;
	SimStates.SGCchangroup = OFF;
	SimStates.SGCchanprams = OFF;
	SimStates.SGCbfh_mode = OFF;
	SimStates.SGCA_mode = OFF;
	SimStates.binary_out = OFF;
	SimStates.gsection = OFF;
	SimStates.binarystartfile = OFF;
	SimStates.startelev = OFF;
	SimStates.latlong = OFF;
	SimStates.dist_routing = OFF;
	SimStates.SGCvoutput = OFF; // switch for sub-grid channel velocity output
	SimStates.DamMode = OFF;
	SimStates.DammaskRead = OFF;
	SimStates.saveint_max = OFF;
	SimStates.maxint = OFF;
	SimStates.ChanMaskRead = OFF;
	SimStates.LinkListRead = OFF;
	SimStates.cuda = OFF;
	SimStates.fv1 = OFF;
	SimStates.fv2 = OFF;
	SimStates.acc_nugrid = OFF;
	SimStates.mwdg2 = OFF;
	SimStates.hwfv1 = OFF;
	SimStates.dg2 = OFF;
	SimStates.dynamicrainfall = OFF;

	SGCchanprams.NSGCprams = 0;
	SGCchanprams.SGCbetahmin = C(0.2);

	/*default resrootname*/
	strcpy(ParFp.res_dirname, "");
	strcpy(ParFp.res_prefix, "res");

	int verbosemode = ReadVerboseMode(argc, argv);

	printversion(verbosemode);

	// if user only wants to know the version then exit
	for (i = 1; i < argc; i++) if (!strcmp(argv[i], "-version")) return(0);

	for (i = 1; i < argc; i++) if (!strcmp(argv[i], "-compare_results"))
	{
		char compare_results_dir1[512];
		char compare_results_dir2[512];
		char compare_results_ext[512];

		if (argc > i + 3)
		{
			sscanf(argv[i + 1], "%511s", compare_results_dir1);
			sscanf(argv[i + 2], "%511s", compare_results_dir2);
			sscanf(argv[i + 3], "%511s", compare_results_ext);
		}
		else
		{
			printf("invalid compare_results options, expect: -compare_results <dirroot> <dirroot> <suffix>\n");
			exit(0);
		}

		for (i = 1; i < argc - 1; i++) if (!strcmp(argv[i], "-resroot"))
		{
			sscanf(argv[i + 1], "%s", ParFp.res_prefix);
			if (verbosemode == ON) printf("Results root name reset by command line: %s\n", ParFp.res_prefix);
		}

		compare_grids(compare_results_dir1, compare_results_dir2, ParFp.res_prefix, compare_results_ext);
		exit(0);
	}

	ReadConfiguration(argc, argv, Fnameptr, Statesptr, Parptr, Solverptr,
			verbosemode);

	// use output folder if requested in parameter file or command line
	if (strlen(ParFp.res_dirname) > 0)
	{
		if (fexist(ParFp.res_dirname) == 0) // check if it doesn't exist
		{
			//create output folder
			sprintf(tmp_sys_com, "%s%s", "mkdir ", ParFp.res_dirname);
			system(tmp_sys_com);
		}
		//set the resroot to include the folder information
		sprintf(ParFp.resrootname, "%s" FILE_SEP"%s", ParFp.res_dirname, ParFp.res_prefix);
	}
	else
	{
		//set to res_prefix
		sprintf(ParFp.resrootname, "%s", ParFp.res_prefix);
	}

	// (MT) redirect all console output to logfile if requested
	if (Statesptr->logfile == ON)
	{
		sprintf(Fnameptr->logfilename, "%s%s", Fnameptr->resrootname, ".log");  //Default log filename
		printf("Redirecting all console output to %s\n\n", Fnameptr->logfilename);
		printf("Lisflood is running ......\n");
		freopen(Fnameptr->logfilename, "w", stdout); // redirect stdout to log file
		setvbuf(stdout, NULL, _IONBF, 0); // set buffer to zero so log file is always up to date
		printversion(verbosemode); // output version here as well (so we get on screen and in file)
	}

	// allow output folder to be determined by commandline
   // for (i = 1; i<argc - 1; i++) if (!strcmp(argv[i], "-dir") || !strcmp(argv[i], "-dirroot"))
   // {
	  //sscanf(argv[i + 1], "%s", ParFp.res_dirname);
   //   SimStates.out_dir=ON;
	  //if (verbosemode == ON) printf("Output folder set by command line: %s\n", ParFp.res_dirname);
   // }
   // // TF: allow results root to be determined by commandline
   // for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-resroot"))
   // {
	  //sscanf(argv[i + 1], "%s", ParFp.res_prefix);
	  //if (verbosemode == ON) printf("Results root name reset by command line: %s\n", ParFp.res_prefix);
   // }

	//// PB: switch to acceleration version
	//for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-acceleration"))
	//{
	   // SimStates.acceleration=ON;
	   // SimStates.adaptive_ts=OFF;
	   // SimStates.qlim=OFF;
	   // if(verbosemode==ON) printf("\nUsing acceleration formulation for floodplain flow\n");
	//}

	//for(i=1;i<argc;i++) if(!strcmp(argv[i],"-checkpoint")) SimStates.checkpoint=ON;

	// A different sim_time if requested
	//for(i=1;i<argc;i++) if(!strcmp(argv[i],"-simtime")) sscanf(argv[i+1],"%" NUM_FMT"",&ParSolver.Sim_Time);

	if (strlen(ParFp.checkpointfilename) == 0)
		sprintf(ParFp.checkpointfilename, "%s%s", ParFp.resrootname, ".chkpnt");  //Default checkpoint filename

	  //for(i=1;i<argc;i++)	if(!strcmp(argv[i],"-loadcheck"))
	  //{
	  //  strcpy(ParFp.loadCheckpointFilename,argv[i+1]);
	  //  SimStates.checkpoint=ON;
	  //  SimStates.checkfile=ON;
	  //}
	if (SimStates.checkpoint == ON && verbosemode == ON)
		printf("Running in checkpointing mode: frequency %" NUM_FMT" hours\n", Params.checkfreq);

	if (SimStates.steadycheck == ON) {
		//for(i=1;i<argc-1;i++) if(!strcmp(argv[i],"-steadytol")) sscanf(argv[i+1],"%" NUM_FMT"",&Params.steadyQtol); // optional tolerance
		if (verbosemode == ON) printf("\nWARNING: simulation will stop on steady-state (tolerance: %.6" NUM_FMT"), or after %.1" NUM_FMT"s.\n", Params.steadyQtol, ParSolver.Sim_Time);
	}

	//code to load in alternative ASCII header for output files
	if (SimStates.alt_ascheader == ON) {
		Params.ascheader = new char*[6];//6 lines in the file
		tmp_fp = fopen(ParFp.ascheaderfilename, "r");
		for (i = 0; i < 6; i++) {
			Params.ascheader[i] = new char[256];//255 characters per line
			fgets(Params.ascheader[i], 255, tmp_fp);
		}
		if (verbosemode == ON) printf("Using alternative ASCII header for output\n");
		fclose(tmp_fp);
	}

	// get system time and echo for user
	if (verbosemode == ON) {
		time_t ts = time(0);
		tm timeS = *localtime(&ts);
		printf("\nStart Date: %d/%d/%d \n", timeS.tm_mday, timeS.tm_mon + 1, timeS.tm_year + 1900);
		printf("Start Time: %d:%d:%d \n\n", timeS.tm_hour, timeS.tm_min, timeS.tm_sec);
	}


	LoadDEM(Fnameptr, Statesptr, Parptr, Arrptr, verbosemode);
	
#ifdef CUDA
	if (Statesptr->cuda == ON)
	{
		if (Statesptr->fv1 == ON)
		{
			lis::cuda::fv1::Simulation simulation;
			simulation.run(ParFp, SimStates, Params, ParSolver, verbosemode);
		}
		else if (Statesptr->fv2 == ON)
		{
			lis::cuda::fv2::Simulation simulation;
			simulation.run(ParFp, SimStates, Params, ParSolver, verbosemode);
		}
		else if (Statesptr->dg2 == ON)
		{
			lis::cuda::dg2::Simulation simulation;
			simulation.run(ParFp, SimStates, Params, ParSolver, verbosemode);
		}
		else if (Statesptr->acceleration == ON)
		{
			lis::cuda::acc::Simulation simulation;
			simulation.run(ParFp, SimStates, Params, ParSolver, verbosemode);
		}
		else if (Statesptr->acc_nugrid == ON)
		{
			lis::cuda::acc_nugrid::Simulation simulation;
			simulation.run(ParFp, SimStates, Params, ParSolver, verbosemode);
		}
		else if (Statesptr->mwdg2 == ON || Statesptr->hwfv1 == ON)
		{
			lis::cuda::adaptive::Simulation simulation;
			simulation.run(argc, argv);
		}
		else
		{
			fprintf(stderr, "cuda only available for acc, fv1, fv2 and dg2 solvers\n");
			return -1;
		}
		return 0;
	}
#endif

	// Dammask needs to be read after LoadDEM and before SGC FEOL
	if (SimStates.DamMode == ON)LoadDamPrams(Fnameptr, Statesptr, Parptr, Damptr, verbosemode); //FEOL
	Damptr->DamLoss = C(0.0); // To ensure dam loss is zero if no dams for mass balance! FEOL
	if (SimStates.DammaskRead == ON)LoadDamMask(Fnameptr, Parptr, Arrptr, Damptr, verbosemode);
	CalcArrayDims(Statesptr, Parptr, Arrptr); // CCS populates dx, dy and dA arrays (calcs correct dimensions if using lat long grid) 

	// dhlin value calculated "on the fly" as a function of dx and gradient (C(0.0002)) from Cunge et al. 1980
	if (SimStates.dhoverw == OFF) ParSolver.dhlin = Params.dx*C(0.0002);

	LoadRiverNetwork(Fnameptr, Statesptr, Parptr, ChannelSegmentsVecPtr, Arrptr, QID7_Vec_Ptr, RiversIndexVecPtr, verbosemode); // CCS
	if (SimStates.ChannelPresent == OFF) ChannelSegments.resize(1); // temp fix to prevent visual studio debuger exiting on the next line (JCN)

	ChannelSegmentType *CSTypePtr = &ChannelSegments[0]; // CCS has to be defined after LoadRiverNetwork has completed.
	int *RiversIndexPtr = &RiversIndex[0];  // CCS has to be defined after LoadRiverNetwork has completed.

	if (QID7.size() != 0) // CCS If there are any tribs then we need to copy the terms from the temp store to the correct place.
	{
		QID7_Store *QID7Ptr = &QID7[0]; // CCS
		UpdateChannelsVector(Statesptr, CSTypePtr, QID7_Vec_Ptr, QID7Ptr, RiversIndexPtr); // CCS
	}

	//override river file friction if specified on command line
	for (i = 1; i < argc - 1; i++) if (!STRCMPi(argv[i], "-nch")) {
		sscanf(argv[i + 1], "%" NUM_FMT"", &tmp);
		if (verbosemode == ON) printf("Channel friction reset by command line: %" NUM_FMT"\n\n", tmp);
		for (chseg = 0; chseg < CSTypePtr->N_Channel_Segments; chseg++) for (i = 0; i < CSTypePtr[chseg].chsz; i++) CSTypePtr[chseg].ChanN[i] = tmp;
	}
	if (SimStates.ChannelPresent == ON) SmoothBanks(Parptr, Solverptr, CSTypePtr, Arrptr, ChannelSegmentsVecPtr, verbosemode);

	if (SimStates.SGC == ON) LoadSGC(Fnameptr, Parptr, Arrptr, Statesptr, verbosemode); // load sub grid channels
	if (SimStates.SGC == ON && SimStates.SGCchanprams == ON) LoadSGCChanPrams(Fnameptr, Statesptr, Parptr, SGCptr, verbosemode); // This loads the parameters for the SGC group information
	if (SimStates.SGC == ON) CalcSGCz(Fnameptr, Statesptr, Parptr, Arrptr, SGCptr, verbosemode);

	if (SimStates.startfile == ON)
    {
        LoadStart(Fnameptr, Statesptr, Parptr, Arrptr, SGCptr, verbosemode);
        if (SimStates.startq2d == ON)
        {
            LoadStartQ2D(Fnameptr, Parptr, Arrptr, verbosemode);
        }
    }
	if (SimStates.binarystartfile == ON) LoadBinaryStart(Fnameptr, Statesptr, Parptr, Arrptr, SGCptr, verbosemode);

	LoadBCs(Fnameptr, Statesptr, Parptr, BCptr, verbosemode);
	LoadBCVar(Fnameptr, Statesptr, Parptr, BCptr, CSTypePtr, Arrptr, ChannelSegmentsVecPtr, verbosemode);
	LoadManningsn(Fnameptr, Parptr, Arrptr, verbosemode);
	LoadDistInfil(Fnameptr, Parptr, Arrptr, verbosemode);
	LoadSGCManningsn(Fnameptr, Parptr, Arrptr, verbosemode);
	// PFU add SGC dirn array
	LoadSGCdirn(Fnameptr, Parptr, Arrptr, verbosemode);
	LoadPor(Fnameptr, Statesptr, Parptr, Arrptr, verbosemode);
	LoadWeir(Fnameptr, Statesptr, Parptr, Arrptr, verbosemode);
	if (SimStates.calc_evap == ON) LoadEvap(Fnameptr, Arrptr, verbosemode);
	if (SimStates.rainfall == ON) LoadRain(Fnameptr, Arrptr, verbosemode);
	if (SimStates.rainfallmask == ON) LoadRainmask(Fnameptr, Parptr, Arrptr, Statesptr, verbosemode);
	if (SimStates.save_stages == ON) LoadStages(Fnameptr, Statesptr, Parptr, Stageptr, verbosemode);
	if (SimStates.gsection == ON) LoadGauges(Fnameptr, Statesptr, Parptr, Stageptr, verbosemode);

	//FEOL note this modifies the DEM! Changes DEM to DEM_NO_DATA where mask is negative
	if (SimStates.routing == ON) // Call FlowDirDEM to generate flow direction map from DEM before main loop CCS
	{
		FlowDirDEM(Parptr, Arrptr, Statesptr, BCptr);
		if (verbosemode == ON) printf("Flow direction map generated from DEM\n\n");
	}

	// apply different starting methods for channel
	if (SimStates.ChannelPresent == ON)
	{
		// calc initial steady state flows down channel
		CalcChannelStartQ(Statesptr, Parptr, Arrptr, CSTypePtr, RiversIndexVecPtr, RiversIndexPtr);

		if (SimStates.startfile == ON)
		{
			// start file is specified. Do nothing, as starting H values for channel already read in from the startfile.
		}
		else if (SimStates.startq == ON)
		{
			// Kinematic: Uses the kinematic initial solution to calculate H from Q
			// Diffusive: Uses diffusive steady state initial solution (default) or can use full dynamic steady state
			// initial if turned on using -dynsw on command line or "ch_dynamic" in the parameter file

			// use the flows to calculate a starting H
			SetChannelStartHfromQ(Statesptr, Parptr, Arrptr, CSTypePtr, Solverptr, RiversIndexVecPtr, RiversIndexPtr);
		}
		else
		{
			// set channel start H to default or user defined H
			SetChannelStartH(Statesptr, Parptr, Arrptr, CSTypePtr, RiversIndexVecPtr, RiversIndexPtr);
		}
	}
	// apply hot starting methods to SGC model
	if (SimStates.startq == ON && SimStates.SGC == ON)
	{
		SGC_hotstart(Statesptr, Parptr, Solverptr, Arrptr);
		if (verbosemode == ON) printf("\nStartq for SGC model implemented\n");
	}

	if (verbosemode == ON) if (SimStates.calc_infiltration == ON) printf("Floodplain infiltration set at: %.10" NUM_FMT" ms-1\n\n", Params.InfilRate);

	//get multiple overpass timings from file
	if (SimStates.multi_op == ON) {
		tmp_fp = fopen(ParFp.opfilename, "r");
		if (tmp_fp != NULL)
		{
			fscanf(tmp_fp, "%i", &Params.op_multinum);
			if (verbosemode == ON) printf("\nMultiple overpass files to be output: %d\n", Params.op_multinum);
			Params.op_multisteps = memory_allocate_numeric_legacy(Params.op_multinum);
			Params.op_multiswitch = new int[Params.op_multinum];
			for (i = 0; i < Params.op_multinum; i++) {
				if (fscanf(tmp_fp, "%" NUM_FMT"", &Params.op_multisteps[i]) != 1) // read in value and check if one value read in successfully
				{
					printf("\nWARNING: overpass file read error at line %i\n", i + 1);
					Params.op_multinum = i; // reset to number of values actually read in
					break;
				}
				Params.op_multiswitch[i] = 0;
				if (verbosemode == ON) printf("Overpass %d at %" NUM_FMT" seconds\n", i, Params.op_multisteps[i]);
			}
			fclose(tmp_fp);
		}
		else {
			SimStates.multi_op = OFF;
			if (verbosemode == ON) printf("\nUnable to open multiple overpass output file: %s\n", ParFp.opfilename);
		}
	}

	//Load checkpointed data if this job has been restarted
	if (SimStates.checkpoint == ON) {
		ReadCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr, CSTypePtr, Arrptr, verbosemode);
		if (verbosemode == ON) printf(" - checkpoint output file: %s\n", ParFp.checkpointfilename);
	}

	//mass balance
	sprintf(t1, "%s%s", ParFp.resrootname, ".mass");
	if (SimStates.checkpoint == ON && ParSolver.t > 0) { //if this is a checkpointed job, we only need to amend the .mass file
		Fps.mass_fp = fopen(t1, "a");
	}
	else {
		Fps.mass_fp = fopen(t1, "w");
	}
	if (Fps.mass_fp != NULL)
	{
		if (ParSolver.t == 0) fprintf(Fps.mass_fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-(Inf+Evap)\n");
		else
		{
			// make a note in the mass file that this is a restart point - user can then edit the overlap out if they want a continuous mass file record.
			fprintf(Fps.mass_fp, "####################################################### Checkpoint restart ########################################################\n");
			fprintf(Fps.mass_fp, "Time         Tstep      MinTstep   NumTsteps    Area         Vol         Qin         Hds        Qout          Qerror       Verror       Rain-(Inf+Evap)\n");
			fflush(Fps.mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
		}
	}
	else
	{
		if (verbosemode == ON)
		{
			printf("Unable to open mass balance file: %s", t1);
			exit(0);
		}
	}
	// FEOL Dam Output file
	//mass balance
	if (SimStates.DamMode == ON)
	{
		sprintf(t1, "%s%s", ParFp.resrootname, ".dam");
		Fps.dam_fp = fopen(t1, "w");

		if (Fps.dam_fp != NULL)
		{
			if (ParSolver.t == 0) fprintf(Fps.dam_fp, "*Time         *Tstep    Area         Vol         Vin         Hds        Vout          Qspill       Qoperation	Rain+Evap)\n");
			else
			{
				fprintf(Fps.dam_fp, "*Time         *Tstep    Area         Vol         Vin         Hds        Vout          Qspill       Qoperation	Rain+Evap)\n");
				fflush(Fps.dam_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
			}
		}
		else
		{
			if (verbosemode == ON)
			{
				printf("Unable to open Dam output file: %s", t1);
				exit(0);
			}
		}
	}

	//stage output file
	if (SimStates.save_stages == ON) {
		sprintf(t1, "%s%s", ParFp.resrootname, ".stage");
		if (SimStates.checkpoint == ON && ParSolver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
			Fps.stage_fp = fopen(t1, "a");
		}
		else {
			Fps.stage_fp = fopen(t1, "w");
		}
		if (Fps.stage_fp != NULL)
		{
			if (ParSolver.t == C(0.0) || SimStates.checkpoint == OFF) // chnage to export z information if initial stage used (JCN)
			{
				fprintf(Fps.stage_fp, "Stage output, depth (m). Stage locations from: %s\n\n", ParFp.stagefilename);
				fprintf(Fps.stage_fp, "Stage information (stage,x,y,elev):\n");
				for (i = 0; i < OutLocs.Nstages; i++)
				{
					if (Statesptr->SGC == ON && Raster.SGCwidth[OutLocs.stage_grid_x[i] + OutLocs.stage_grid_y[i] * Parptr->xsz] > 0) // if a SUB GRID channel is present export the channel bed elevation
					{
						if (OutLocs.stage_check[i] == 1) fprintf(Fps.stage_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\t%.4" NUM_FMT"\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i], Raster.SGCz[OutLocs.stage_grid_x[i] + OutLocs.stage_grid_y[i] * Parptr->xsz]);
						else fprintf(Fps.stage_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\tn/a\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i]);
					}
					else
					{
						if (OutLocs.stage_check[i] == 1) fprintf(Fps.stage_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\t%.4" NUM_FMT"\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i], Raster.DEM[OutLocs.stage_grid_x[i] + OutLocs.stage_grid_y[i] * Parptr->xsz]);
						else fprintf(Fps.stage_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\tn/a\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i]);
					}
				}
				fprintf(Fps.stage_fp, "\nOutput, depths:\n");
				fprintf(Fps.stage_fp, "Time; stages 1 to %d\n", OutLocs.Nstages);
			}
			else
			{
				fprintf(Fps.stage_fp, "####################################################### Checkpoint restart ########################################################\n");
				fflush(Fps.stage_fp);
			}
		}
		else
		{
			if (verbosemode == ON) printf("Unable to open stage output file: %s", t1);
			SimStates.save_stages = OFF;
		}

	}
	//velocity output file
	if (SimStates.save_stages == ON && Statesptr->voutput_stage == ON)
	{
		sprintf(t1, "%s%s", ParFp.resrootname, ".velocity");
		if (SimStates.checkpoint == ON && ParSolver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
			Fps.vel_fp = fopen(t1, "a");
		}
		else {
			Fps.vel_fp = fopen(t1, "w");
		}
		if (Fps.vel_fp != NULL) {
			if (ParSolver.t == 0) {
				fprintf(Fps.vel_fp, "Velocity output, velocity (ms-1). Velocity locations from: %s\n\n", ParFp.stagefilename);
				fprintf(Fps.vel_fp, "Stage information (stage,x,y,elev):\n");
				for (i = 0; i < OutLocs.Nstages; i++) {
					if (OutLocs.stage_check[i] == 1) fprintf(Fps.vel_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\t%.4" NUM_FMT"\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i], Raster.DEM[OutLocs.stage_grid_x[i] + OutLocs.stage_grid_y[i] * Params.xsz]);
					else fprintf(Fps.vel_fp, "%d\t%.4" NUM_FMT"\t%.4" NUM_FMT"\tn/a\n", i + 1, OutLocs.stage_loc_x[i], OutLocs.stage_loc_y[i]);
				}
				fprintf(Fps.vel_fp, "\nOutput, depths:\n");
				fprintf(Fps.vel_fp, "Time; velocities 1 to %d\n", OutLocs.Nstages);
			}
			else {
				fprintf(Fps.vel_fp, "####################################################### Checkpoint restart ########################################################\n");
				fflush(Fps.vel_fp);
			}
		}
		else {
			if (verbosemode == ON) printf("Unable to open velocity output file: %s", t1);
			SimStates.save_stages = OFF;
		}
	}

	//velocity output file
	if (SimStates.gsection == ON)
	{
		sprintf(t1, "%s%s", ParFp.resrootname, ".discharge");
		if (SimStates.checkpoint == ON && ParSolver.t > 0) { //if this is a checkpointed job, we only need to amend the .stage file
			Fps.gau_fp = fopen(t1, "a");
		}
		else {
			Fps.gau_fp = fopen(t1, "w");
		}
		if (Fps.gau_fp != NULL) {
			if (ParSolver.t == 0) {
				fprintf(Fps.gau_fp, "Discharge output, discharge (m3s-1). Discharge locations from: %s\n\n", ParFp.gaugefilename);
				fprintf(Fps.gau_fp, "Time; discharge 1 to %d\n", OutLocs.Ngauges);
			}
			else {
				fprintf(Fps.gau_fp, "####################################################### Checkpoint restart ########################################################\n");
				fflush(Fps.gau_fp);
			}
		}
		else {
			if (verbosemode == ON) printf("Unable to open discharge output file: %s", t1);
			SimStates.gsection = OFF;
		}
	}

	////find out if we are going to compress output on the fly
	//for(i=1;i<argc;i++) {
	//  if(!strcmp(argv[i],"-gzip")) {
	//    SimStates.call_gzip=ON;
	   // SimStates.output_params.call_gzip = ON;
	//    if(verbosemode==ON) printf("\nOutput will be compressed using Gzip\n");
	//  }
	//}
	if (Statesptr->maxdepthonly == ON)
	{
	}
	else
	{
		// output debug files (used DEM, channel mask seg mask) if required
		if (Statesptr->debugmode == ON)
			debugfileoutput(Fnameptr, Statesptr, Parptr, Arrptr);

		if (Statesptr->SGC == ON) // output base/bed DEM including channel depths for display purposes with water depth
			write_ascfile(Fnameptr->resrootname, -1, ".dem", Arrptr->SGCz, Arrptr->DEM, 0, Statesptr, Parptr);
		else  // Write out final DEM if not subgrid - includes 1D river channel and channel bank modifications
			write_ascfile(Fnameptr->resrootname, -1, ".dem", Arrptr->DEM, Arrptr->DEM, 0, Statesptr, Parptr);
	}

	//start simulation
	time(&ParSolver.time_start);
	if (Statesptr->SGC == ON) // SGC output
	{
		Fast_MainStart(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr, Stageptr, CSTypePtr, Arrptr, SGCptr, ChannelSegmentsVecPtr, Damptr, verbosemode); //Damptr added by FEOL
		//IterateQ(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr, Stageptr, CSTypePtr, Arrptr, SGCptr, RiversIndexVecPtr, RiversIndexPtr, ChannelSegmentsVecPtr, verbose);
	}
	else if (Statesptr->fv1 == ON)
	{
		fv1::solve(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr,
				Stageptr, Arrptr, verbosemode);
	}
	else if (Statesptr->dg2 == ON)
	{
		dg2::solve(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr,
				Stageptr, Arrptr, verbosemode);
        //dg2new::DG2Solver solver(Fnameptr, &Fps, Statesptr, Parptr, Solverptr,
        //        BCptr, Stageptr, Arrptr, verbosemode);
        //solver.solve();

        /*
		dg2::solve(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr,
				Stageptr, Arrptr, verbosemode);
        */
	}
	else
	{
		IterateQ(Fnameptr, &Fps, Statesptr, Parptr, Solverptr, BCptr, Stageptr, CSTypePtr, Arrptr, SGCptr, RiversIndexVecPtr, RiversIndexPtr, ChannelSegmentsVecPtr, verbosemode);
	}
	time(&ParSolver.time_finish);

	//Final checkpoint
	if (SimStates.checkpoint == ON) WriteCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr, CSTypePtr, Arrptr, verbosemode);

	// get system time and echo for user
	if (verbosemode == ON) {
		time_t tf = time(0);
		tm timeF = *localtime(&tf);
		printf("\nFinish Date: %d/%d/%d \n", timeF.tm_mday, timeF.tm_mon + 1, timeF.tm_year + 1900);
		printf("Finish Time: %d:%d:%d \n\n", timeF.tm_hour, timeF.tm_min, timeF.tm_sec);
	}

	//iteration time
	ParSolver.itrn_time = ParSolver.itrn_time + (NUMERIC_TYPE)difftime(ParSolver.time_finish, ParSolver.time_start);
	if (verbosemode == ON) printf("\n  Total computation time: %.2" NUM_FMT" mins\n\n", (ParSolver.itrn_time / C(60.0)));

	if (SimStates.logfile == ON)
	{
		freopen("CON", "w", stdout);
		printf("\nLisflood run finished see log file for run details");
	}

	if (SimStates.save_stages == ON) fclose(Fps.stage_fp);
	sprintf(t1, "%s%s", ParFp.resrootname, ".stage");
	if (SimStates.call_gzip == ON) {
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", t1);
		system(tmp_sys_com);
	}

	if (Fps.mass_fp != NULL) fclose(Fps.mass_fp);
	sprintf(t1, "%s%s", ParFp.resrootname, ".mass");
	if (SimStates.call_gzip == ON) {
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", t1);
		system(tmp_sys_com);
	}
	// FEOL Dam output
	if (Statesptr->DamMode == ON)
	{
	
	if (Fps.dam_fp != NULL) fclose(Fps.dam_fp);
	sprintf(t1, "%s%s", ParFp.resrootname, ".dam");
	}
  
  
  return 0;
}
//---------------------------------------------------------------------------
void printversion(int verbose)
// printout header with program and version number
{
  printf("***************************\n");
  printf(" LISFLOOD-FP version %d.%d.%d (%s)\n", LF_VersionMajor, LF_VersionMinor, LF_VersionInc, NUMERIC_TYPE_NAME);
  if (verbose == ON)
  {
#if defined (__INTEL_COMPILER)
	  printf("Intel Compiler version: %d\n", __INTEL_COMPILER);

	  //https://software.intel.com/en-us/node/514528
	  printf("CPU instructions used:");
#if defined (__AVX2__)
	  printf(" AVX2");
#endif
#if defined (__AVX__)
	  printf(" AVX");
#endif
#if defined (__SSE4_2__)
	  printf(" SSE_4.2");
#endif
#if defined (__SSE4_1__)
	  printf(" SSE_4.1");
#endif
#if defined (__SSE3__)
	  printf(" SSE3");
#endif
#if defined (__SSE2__)
	  printf(" SSE2");
#endif
#if defined (__SSE__)
	  printf(" SSE");
#endif
	  printf("\n");

#endif
  }

#if defined (CUDA)
  printf("CUDA supported\n");
#endif
#if defined (_PROFILE_MODE) && _PROFILE_MODE > 0
  printf("Profile Mode Enabled: %d\n", _PROFILE_MODE);
#endif
#if defined (_SGM_BY_BLOCKS) && _SGM_BY_BLOCKS > 0
  printf("_SGM_BY_BLOCKS: %d\n", _SGM_BY_BLOCKS);
#endif
#if defined (_BALANCE_TYPE) && _BALANCE_TYPE > 0
  printf("_BALANCE_TYPE: %d\n", _BALANCE_TYPE);
#endif
#if defined (_ONLY_RECT) && _ONLY_RECT == 1
  printf("Rectangular channels only.\n");
#endif
#if defined (_DISABLE_WET_DRY) && _DISABLE_WET_DRY == 1
  printf("_DISABLE_WET_DRY.\n");
#endif
#if defined (_CALCULATE_Q_MODE) && (_CALCULATE_Q_MODE != 0)
  printf("_CALCULATE_Q_MODE %d.\n", _CALCULATE_Q_MODE);
#endif

  printf("***************************\n\n");
}

