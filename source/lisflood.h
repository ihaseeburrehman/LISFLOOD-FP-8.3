#pragma once
/*
#####################################################################################
LISFLOOD-FP flood inundation model
#####################################################################################

� copyright Bristol University Hydrology Research Group 2008

webpage -	http://www.ggy.bris.ac.uk/research/hydrology/models/lisflood
contact -	Professor Paul Bates, email: paul.bates@Bristol.ac.uk,
Tel: +44-117-928-9108, Fax: +44-117-928-7878

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <omp.h>
#include <vector> // CCS
#include <iostream> // CCS
//#include <netcdf.h> // JCN

#ifndef _SGM_BY_BLOCKS
// 0 SGC by row 
// 1 SGC by block
#define _SGM_BY_BLOCKS 0
#endif
#ifndef _BALANCE_TYPE
// 1 balance by wet, 
// 0 even fixed balance
#define _BALANCE_TYPE 0 
#endif

#ifndef _NUMERIC_MODE
#define _NUMERIC_MODE 1
#endif
#ifndef _PROFILE_MODE
#define _PROFILE_MODE 0
#endif
#ifndef _ONLY_RECT
#define _ONLY_RECT 1
#endif
#ifndef _DISABLE_WET_DRY
#define _DISABLE_WET_DRY 0
#endif
#ifndef _CALCULATE_Q_MODE
#define _CALCULATE_Q_MODE 1
#endif
#ifndef _NETCDF
#define _NETCDF 1
#endif

//#define TESTING
#ifdef TESTING
void RunTests();
#endif
#ifdef _MSC_VER
//#define RESULT_CHECK 1
#endif

#ifdef _DEBUG
//#define RESULT_CHECK 1
#endif

// older versions of visual studio does not contain cbrt function
#ifndef cbrt
	#define cbrt(x) pow(x,1/3.0);
#endif
#ifndef cbrtf
	#define cbrtf(x) powf(x,1/3.0f);
#endif

#if _NUMERIC_MODE == 1
#define NUMERIC_TYPE double
#define NUMERIC_TYPE_NAME "double"
#define NUM_FMT "lf"
#define FMAX fmax
#define FMIN fmin
#define POW pow
#define FABS fabs
#define CBRT cbrt
#define SQRT sqrt
#define C(x) x
#else
#define NUMERIC_TYPE float
#define NUMERIC_TYPE_NAME "float"
#define NUM_FMT "f"
#define FMAX fmaxf
#define FMIN fminf
#define POW powf
#define FABS fabsf
#define CBRT cbrtf
#define SQRT sqrtf
#define C(x) x##f
#endif

#ifdef __unix__
#define FILE_SEP "/"
#define STRCMPi strcasecmp
#elif __APPLE__
#define FILE_SEP "/"
#define STRCMPi strcasecmp
#else
#define STRCMPi strcmpi
#define FILE_SEP "\\"
#endif

#if _XOPEN_SOURCE >= 600 || _ISOC99_SOURCE || _POSIX_C_SOURCE >= 200112L || _MSC_VER >= 1800
#define getmax(a, b) FMAX(a, b)
#else
#define getmax(a, b) (a>b?a:b)
#endif

#if _XOPEN_SOURCE >= 600 || _ISOC99_SOURCE || _POSIX_C_SOURCE >= 200112L || _MSC_VER >= 1800
#define getmin(a, b) FMIN(a, b)
#else
#define getmin(a, b) (a<b?a:b)
#endif

#define LINE_BUFFER_LEN 4096

// define basic constants
#define ON 1
#define OFF 0
#define CHKINTERVAL C(1.0) // default checkpoint interval in hours runtime
#define NULLVAL C(-9999.0) // MT: define ascii file NULL value as constant
#define DEM_NO_DATA C(1e7) // new code should use Pars.nodata_elevation instead // decreased to 1e7 to support float
#define DEFAULT_PRECISION 6

#define PARAM_FILE 0
#define CMD_LINE 1

#define TIME 1
#define TIME_SPACE 2

#define fix_small_negative(value) ((value < C(0.0) && value > -0.0001) ? C(0.0) : value)

#pragma warning( disable : 4996)  // MT: disable visual C++ depreciation warnings so we can see other warnings
// #pragma warning( disable : 1478)  // MT: disable intel depreciation warnings so we can see other warnings [CCS VS2010 throws warning C4616: #pragma warning : warning number '1478' not a valid compiler warning]

using namespace std; // CCS

/*
*****************************************************************************

Define the structures
---------------------
Fnames - Contains all the filenames from the .par file
States - Contains all the state parameters for the simulation
Pars - Contains the parameter values specified in the .par file
Solver - Defines the solution settings
BoundaryValues - Used in the DG2 SWE solver
Arrays - Defines the global arrays
ChannelSegmentType - Defines the channel
Stage - Variables for outputting stage information at specified locations
Files - General output file pointers
BoundCs - Boundary conditions

*****************************************************************************
*/

/// time series loaded from .bdy file
struct TimeSeries{
	NUMERIC_TYPE *time;
	NUMERIC_TYPE *value;
	int count;
	int prev_index;

	// store the prev time queried
	NUMERIC_TYPE prev_time;
	// store the prev time queried result value
	NUMERIC_TYPE prev_value;
};

enum ESourceType {
	NONE0 = 0,
	FREE1 = 1,
	HFIX2 = 2,
	HVAR3 = 3,
	QFIX4 = 4,
	QVAR5 = 5,
	// FREE or Qout(rivers)
	FREE6 = 6,
	// rivers
	TRIB7 = 7,
	// rivers
	RATE8 = 8,
};

enum EWeirType
{
	EWeir_Weir = 0,
	EWeir_Bridge = 1
};

enum EDirection
{
	//N = 1, E = 2, S = 3, W = 4
	DirectionNA = 0,
	North = 1,
	East = 2,
	South = 3,
	West = 4
};

typedef struct
{
	NUMERIC_TYPE *DEM;
	NUMERIC_TYPE *H;
	NUMERIC_TYPE *HU;
	NUMERIC_TYPE *HV;
} BoundaryValues;

/*! \struct Arrays
Stores the pointers to arrays required globally in the computation. Defined as 1D
vectors but stores 2D data determined by the array subscripts.
*/
struct Arrays{
	/*! DEM, Water height, Flow in x-direction and Flow in y-direction */
	NUMERIC_TYPE *DEM; // Digital elevation model
	NUMERIC_TYPE *H;
	NUMERIC_TYPE *Qx;
	NUMERIC_TYPE *Qy;
	NUMERIC_TYPE *Qxold;
	NUMERIC_TYPE *Qyold;
	NUMERIC_TYPE *U;
	NUMERIC_TYPE *Rainmask;  //Distrubted rainfall AS
	
	/* Fields specific to Roe and SWE solvers  */
	NUMERIC_TYPE *HU;
	NUMERIC_TYPE *HV;
	NUMERIC_TYPE *FHx;
	NUMERIC_TYPE *FHUx;
	NUMERIC_TYPE *FHVx;
	NUMERIC_TYPE *FHy;
	NUMERIC_TYPE *FHUy;
	NUMERIC_TYPE *FHVy;

	/* Roe-specific fields */
	NUMERIC_TYPE *RSHU;
	NUMERIC_TYPE *LSHU;
	NUMERIC_TYPE *RSHV;
	NUMERIC_TYPE *LSHV;
	NUMERIC_TYPE *BSHU;
	NUMERIC_TYPE *TSHU;
	NUMERIC_TYPE *BSHV;
	NUMERIC_TYPE *TSHV;

	/* SWE-specific fields */
	NUMERIC_TYPE *Zstar_x;
	NUMERIC_TYPE *Zstar_y;
	NUMERIC_TYPE *Hstar_neg_x;
	NUMERIC_TYPE *Hstar_pos_x;
	NUMERIC_TYPE *Hstar_neg_y;
	NUMERIC_TYPE *Hstar_pos_y;

	/* DG2-specific fields */
	NUMERIC_TYPE *DEM1x;
	NUMERIC_TYPE *DEM1y;
	NUMERIC_TYPE *H1x;
	NUMERIC_TYPE *H1y;
	NUMERIC_TYPE *HU1x;
	NUMERIC_TYPE *HU1y;
	NUMERIC_TYPE *HV1x;
	NUMERIC_TYPE *HV1y;

	NUMERIC_TYPE *H_int;
	NUMERIC_TYPE *H1x_int;
	NUMERIC_TYPE *H1y_int;
	NUMERIC_TYPE *HU_int;
	NUMERIC_TYPE *HU1x_int;
	NUMERIC_TYPE *HU1y_int;
	NUMERIC_TYPE *HV_int;
	NUMERIC_TYPE *HV1x_int;
	NUMERIC_TYPE *HV1y_int;

	NUMERIC_TYPE *ETA1x_slopelim;
	NUMERIC_TYPE *ETA1y_slopelim;
	NUMERIC_TYPE *HU1x_slopelim;
	NUMERIC_TYPE *HU1y_slopelim;
	NUMERIC_TYPE *HV1x_slopelim;
	NUMERIC_TYPE *HV1y_slopelim;

	NUMERIC_TYPE *HUstar_neg_x;
	NUMERIC_TYPE *HUstar_pos_x;
	NUMERIC_TYPE *HUstar_neg_y;
	NUMERIC_TYPE *HUstar_pos_y;
	NUMERIC_TYPE *HVstar_neg_x;
	NUMERIC_TYPE *HVstar_pos_x;
	NUMERIC_TYPE *HVstar_neg_y;
	NUMERIC_TYPE *HVstar_pos_y;

	/* Flow direction map for Rainfall*/
	int *FlowDir; // CCS: added to hold DEM flow direction map for routing shallow rainfall flow 13/03/2012
	NUMERIC_TYPE *Route_dH;
	NUMERIC_TYPE *RouteInt; // CCS: added to record routing scheme dH and interval
	

	/* ---------------- */
	NUMERIC_TYPE *maxH;
	NUMERIC_TYPE *maxHtm;
	NUMERIC_TYPE *initHtm;
	NUMERIC_TYPE *totalHtm;
	NUMERIC_TYPE *Manningsn;
	NUMERIC_TYPE *SGCManningsn;
	NUMERIC_TYPE *paerial;
	NUMERIC_TYPE *pbound;

	int weir_count;
	//lists of weir data (each have 'weir_count' items)
	NUMERIC_TYPE *Weir_hc;
	NUMERIC_TYPE *Weir_Cd;
	NUMERIC_TYPE *Weir_m;
	NUMERIC_TYPE *Weir_w;
	EDirection *Weir_Fixdir;
	EWeirType *Weir_Typ;

	//grid for weirs updating Qx
	int *Weir_Identx;
	//grid for weirs updating Qy
	int *Weir_Identy;

	TimeSeries * evap;
	TimeSeries * rain;

	int *ChanMask;
	int *SegMask;

	NUMERIC_TYPE *TRecx;
	NUMERIC_TYPE *TRecy; // MT: add to record TStep

	NUMERIC_TYPE *LimQx;
	NUMERIC_TYPE *LimQy; // MT: add to record Qlimits
	NUMERIC_TYPE *Vx;
	NUMERIC_TYPE *Vy;
	NUMERIC_TYPE *maxVx;
	NUMERIC_TYPE *maxVy;
	NUMERIC_TYPE *Vc; // JCN: added to record velocity
	NUMERIC_TYPE *maxVc;
	NUMERIC_TYPE *maxVcH;
	NUMERIC_TYPE *maxHaz; // JCN added to calculate hazard
	NUMERIC_TYPE *SGCwidth;
	NUMERIC_TYPE *SGCz;
	NUMERIC_TYPE *QxSGold;
	NUMERIC_TYPE *QySGold;
	NUMERIC_TYPE *SGCbfH;
	NUMERIC_TYPE *SGCVol;
	NUMERIC_TYPE *SGCdVol;
	NUMERIC_TYPE *SGCbfV;
	NUMERIC_TYPE *SGCc;
	NUMERIC_TYPE *SGCFlowWidth;
	NUMERIC_TYPE *SGCdx;
	NUMERIC_TYPE *SGCcat_area;// JCN added to store widths and depths
	NUMERIC_TYPE *dx;
	NUMERIC_TYPE *dy;
	NUMERIC_TYPE *dA; // CCS added for lat long data
	NUMERIC_TYPE *DamMask; // FEOL for Res..
	NUMERIC_TYPE *dist_infiltration; // JCN stores distributed infiltration rates
	
	int  *SGCgroup;
	int  *SGCdirn;  // PFU for prescribing sub grid channel flow directions
	BoundaryValues boundary;
};

//-------------------------------------------
// Files
struct Files{
	FILE *mass_fp;
	FILE *stage_fp;
	FILE *vel_fp;
	FILE *gau_fp;
	FILE *dam_fp;
};

//-------------------------------------------
// Fnames
struct Fnames{

	char resrootname[512]; // resrootname will be res_dirname + res_prefix
	char demfilename[256];
	char startfilename[256];
	char chanfilename[256];

	char res_prefix[256];
	char res_dirname[256];
	char qfilename[256];
	char nfilename[256];
	char SGCnfilename[256];
	char porfilename[256];
	char rivername[256];
	char bcifilename[256];
	char bdyfilename[256];
	char weirfilename[256];
	char opfilename[256];
	char stagefilename[256];
	char ascheaderfilename[256];
	char multiriverfilename[256]; // CCS
	char checkpointfilename[256]; // used to write a checkpoint file (note it is placed in the input file dir not the results dir)
	char loadCheckpointFilename[256]; // explicit checkpoint file to start run (specify with -loadcheck option, defaults to checkpointfilename)
	char evapfilename[256];
	char rainfilename[256];
	char rainmaskname[256]; //Distributed rainfall AS
	char logfilename[256];
	char SGCwidthfilename[256]; // JN sub grid channel widths
	char SGCbankfilename[256]; // JN sub grid channel bank elevations
	char SGCbedfilename[256];  // JN sub grid channel bed elevation
	char SGCleveefilename[256];  // NQ levee addition
	char SGCcat_areafilename[256]; // JN sub grid channel accumulation area
	char SGCchangroupfilename[256];
	char SGCchanpramsfilename[256];
	char gaugefilename[256];
	char DamMaskfilename[256]; // FEOL
	char Damparfilename[256]; // FEOL
	char ChanMaskfilename[256]; // JCN
	char LinkListfilename[256]; // JCN
	char SGCdirnfilename[256];  // PFU
	char infilfilename[256]; // JCN
    char dynamicrainfilename[256];
};


//-------------------------------------------
// Boundary Conditions
struct BoundCs{
	int* xpi; //used in legacy and read in
	int* ypi; //used in legacy and read in

	char  *PS_Name;
	ESourceType   *PS_Ident;
	int   numPS;
	// PS_Val used in case of fixed e.g. HFIX or QFIX (otherwise set to -1)
	NUMERIC_TYPE *PS_Val;
	// time series indexed by psi (point source index) //TFD
	// PS_TimeSeries used in case of var e.e. HVAR or QVAR (otherwise set to NULL)
	TimeSeries **PS_TimeSeries;

	NUMERIC_TYPE *PS_Q_FP_old;
	NUMERIC_TYPE *PS_Q_SG_old;

	ESourceType   *BC_Ident;
	int numBCs;
	char  *BC_Name;
	// BC_Val used in case of fixed e.g. HFIX or QFIX (otherwise set to -1)
	NUMERIC_TYPE *BC_Val;
	// time series indexed by bci (boundary condition index) //TFD
	// BC_TimeSeries used in case of var e.e. HVAR or QVAR (otherwise set to NULL)
	TimeSeries **BC_TimeSeries;

	NUMERIC_TYPE Qpoint_pos; // replace Qpoint with positive and negative versions to keep track of input or output for point sources
	NUMERIC_TYPE Qpoint_neg;
	NUMERIC_TYPE Qin;
	NUMERIC_TYPE Qout;
	NUMERIC_TYPE QChanOut;
	NUMERIC_TYPE VolInMT; // added by JCN stores volume in over mass inteval
	NUMERIC_TYPE VolOutMT; // added by JCN stores volume out over mass inteval

	std::vector<TimeSeries> allTimeSeries;
};

//-------------------------------------------
// Stage
struct Stage{
	int Nstages, Ngauges;
	NUMERIC_TYPE *stage_loc_x, *stage_loc_y;
	NUMERIC_TYPE *gauge_loc_x, *gauge_loc_y, *gauge_dist;
	int *stage_grid_x, *stage_grid_y, *stage_check;
	int *gauge_grid_x, *gauge_grid_y;
	int *gauge_cells;
	EDirection *gauge_dir;
};

// SGC parameters
struct SGCprams{
	int NSGCprams;
	int *SGCchantype;
	NUMERIC_TYPE SGCbetahmin;
	NUMERIC_TYPE *SGCp, *SGCr, *SGCs, *SGCm, *SGCa;
	//mannings squared, indexed by channel group
	NUMERIC_TYPE *SGCn;
	NUMERIC_TYPE *SGCgamma, *SGCbeta1, *SGCbeta2, *SGCbeta3, *SGCbeta4, *SGCbeta5;
};

struct NetCDFVariable
{
  int ncid;
  size_t xlen;
  size_t ylen;
  size_t tlen;
  int varid;
  size_t time_idx;
  NUMERIC_TYPE dt;
  NUMERIC_TYPE* times;
  NUMERIC_TYPE* xs;
  NUMERIC_TYPE* ys;
  NUMERIC_TYPE* data;
};

struct NetCDFState
{
	int init_done;

	int ncid;
	int dimid_time;
	int dimid_x;
	int dimid_y;

	int dimid_x_edge; //x+1
	int dimid_y_edge; //y+1

	//int dimid_height;
	//int dimid_q;
	//int dim_id_Velocity;

	// one dimentional - record time at each write
	int varid_time;
	int varid_x;
	int varid_y;

	// time series variables
	int varid_depth;
	int varid_elevation;
	
	int varid_qx;
	int varid_qy;
	int varid_qcx;
	int varid_qcy;
	int varid_Vx;
	int varid_Vy;
	//int varid_Vc; // previously not written
	int varid_sgc_Vx;
	int varid_sgc_Vy;
	int varid_sgc_Vc;


	// single grid variables
	int varid_inittm;
	int varid_totaltm;
	int varid_max;
	int varid_mxe;
	int varid_maxtm;
	int varid_maxVx;
	int varid_maxVy;
	int varid_maxVc;
	int varid_maxVcd;
	int varid_maxHaz;

	

};


struct OutputParams
{
	int standard_extensions;
	int ascii_out;
	int binary_out;

	int call_gzip;

	int netcdf_out;
	NetCDFState netcdf_state;
};


//-------------------------------------------
// Simulation States
struct States{
	int ChannelPresent;
	int TribsPresent;
	int NCFS;
	int save_depth;
	int save_elev;
	int save_vtk;
	int single_op;
	int multi_op;
	int calc_area;
	int calc_meandepth;
	int calc_volume;
	int save_stages;
	int adaptive_ts;
	int acceleration; // PB: Flag to switch to acceleration formulation
	int qlim; // TJF: Flag for qlim version
	int debugmode;
	int save_Qs;
	int calc_infiltration;
	int calc_distributed_infiltration; // uses a file to use spatially distributed infiltration rates
	int call_gzip;
	int alt_ascheader;
	int checkpoint;
	int checkfile;
	int calc_evap;
	int rainfall; // TJF: added for time varying, spatially uniform rainfall
	int rainfallmask; //added for distributed rainfall
	int routing; // CCS: added for routing routine. 
	int routing_mass_check; // CCS: added for routing routine. 
	int diffusive_switch; //TFD switch for routing
	int reset_timeinit;
	int profileoutput;
	int porosity;
	int weirs;
	int save_Ts;   // MT: added flag to output adaptive timestep
	int save_QLs;  // MT: added flag to output Qlimits
	int diffusive; // MT: added flag to indicate wish to use diffusive channel solver instead of default kinematic
	int startq;    // MT: added flag to indicate wish to use start flow to calculate initial water depths throughout channel
	int logfile;   // MT: added flag to record logfile
	int startfile; // MT: added flag to note use of startfile
	int start_ch_h; // MT: added flag to note use of starting H in channel
	int comp_out; // TJF: added to make computational output information optional
	int chainagecalc; // MT: added so user can switch off grid independent river chainage calculation
	int mint_hk; // JN: added to request maxH, maxHtm totalHtm and initHtm be calulated at the mass interval
	int Roe; // JN/IV: added to use Roe solver
	int killsim; // MDW: added to flag kill of simulation after specified run time
	int dhoverw; // TJF: added as a switch for dhlin (ON - dhlin set by command line/parfile; OFF - dhlin prescribed by gradient C(0.0002) Cunge et al. 1980)
	int drychecking; //JN Option to turn DryCheck off
	int voutput; // exports velocity esimates based on Q's of Roe velocity (JCN)
	int maxdepthonly; // only export maxdepth file (AS)
	int voutput_max; // if max not required, v doesn't need to be calculated each time step
	int voutput_stage; // only if stages already enabled - can save velocity with stage, without saving velocity grids
	int steadycheck; // MDW: added flag to check for model steady state
	int hazard; // JN additional module for calculating hazards
	int startq2d; // JN: initalises inertial model with uniform flow for Qold
	int Roe_slow; // JN: ghost cell version of Roe solver
	int multiplerivers; // CCS multiple river switch
	int SGC; // JN sub gird channels r
	int SGCbed; // JN sub grid bed elevation file to override hydraulic geometry
	int SGClevee; // NQ levee (outflow only) addition
	int SGCcat_area; // JN sub grid channel accumulated area to override hydraulic geometry based on width
	int SGCchangroup; // turns on distributed channe groups
	int SGCchanprams; // parameters for distributed channel groups
	int binary_out; // JN binary raster output
	OutputParams output_params;
	int gsection; // JN virtual gauge sections
	int binarystartfile; // JN load a binary start file
	int startelev; // used to use an elevation file for the startfile
	int latlong; // CCS: added for lat-long coordinate systems
	int SGCbfh_mode; // JCN switches model to use parameter p as bank full depth
	int SGCA_mode; // JCN switches model to use parameter p as bank full Area
	int dist_routing; // JCN turnes on spatially distributed routing velocity
	int SGCvoutput; // JCN Turns on sub-grid channel velocity output
	int DamMode; // FEOL Turns on reservoir/dam
	int DammaskRead; // FEOL Turns on reservoir/dam
	int ChanMaskRead; // JCN Read channel mask
	int LinkListRead; // JCN Read channel mask
	int saveint_max; //instructs model to save max depth at each stage interval
	int maxint;
	int SGCd8; //PFU flag to choose d8 directions instead of d4 in the sub grid channels
	int cuda;
	int fv1;
	int fv2;
	int dg2;
	int dynamicrainfall;
	int acc_nugrid;
	int mwdg2;
	int hwfv1;
};


//-------------------------------------------
// Model Parameters
struct Pars{
	int xsz, ysz;
	NUMERIC_TYPE dx, dx_sqrt;
	NUMERIC_TYPE dy, dA;
	// friction flood plain (when per cell mannings disabled) - (not squared)
	NUMERIC_TYPE FPn;
	NUMERIC_TYPE tlx, tly, blx, bly;
	NUMERIC_TYPE SaveInt, MassInt;
	NUMERIC_TYPE SaveTotal, MassTotal;
	int SaveNo;
	int op_multinum;
	NUMERIC_TYPE *op_multisteps;
	int *op_multiswitch;
	NUMERIC_TYPE op;
	NUMERIC_TYPE InfilRate;
	NUMERIC_TYPE InfilLoss, EvapLoss, RainLoss; // previous mass interval loss
	NUMERIC_TYPE InfilTotalLoss, EvapTotalLoss, RainTotalLoss; // cumulative loss
	NUMERIC_TYPE checkfreq, nextcheck;
	NUMERIC_TYPE reset_timeinit_time;
	NUMERIC_TYPE maxelev, zlev; // Water depth dependent porosity
	int zsz; // Water depth dependent porosity
	int Por_Ident;
	NUMERIC_TYPE dAPor;
	char **ascheader;
	NUMERIC_TYPE ch_start_h; // starting depth of channel flow. default to 2m or read from par file.
	NUMERIC_TYPE killsim_time; // time to kill simulation
	NUMERIC_TYPE steadyQdiff, steadyQtol, steadyInt, steadyTotal; // used for checking steady-state
	NUMERIC_TYPE SGC_p; // sub grid channel width depth exponent
	NUMERIC_TYPE SGC_r; // sub grid channel width depth mutiplier
	int SGCchan_type; // JCN bank slop for trapazoidal channel
	NUMERIC_TYPE SGC_s, SGC_2, SGC_n; // JCN trapazodal channel slope dependent constant
	NUMERIC_TYPE *SGCprams; // pointer to table of SGC parameters
	NUMERIC_TYPE Routing_Speed, RouteInt; // CCS variables controlling routing speed in rainfall routing routine
	NUMERIC_TYPE RouteSfThresh; // CCS water surface slope at which routing scheme takes over from shallow water eqn is SGC mode (when routing==ON).
	NUMERIC_TYPE DiffusiveFroudeThresh;
	NUMERIC_TYPE SGC_m, SGC_a; // allows a meander coefficient to be set for the sub-grid model, default 1, allows channel upstream area to be set, defaul -1;
	NUMERIC_TYPE min_dx, min_dy, min_dx_dy; // CCS added to hold minimum values of dx and dy when using lat-long projected grids.
	NUMERIC_TYPE max_Froude; // maximum Froude for sub grid solver (needs CALCULATE_Q_MODE = 1)
	NUMERIC_TYPE maxint; // writes and resets maximum depth over interval 
	NUMERIC_TYPE maxintTotal; // writes and resets maximum depth over interval 
	int maxintcount; // counts number of maxint saves
	int output_precision;
    NUMERIC_TYPE nodata_elevation; // DEM elevation used for NODATA values
    int drain_nodata; // remove water from DEM NODATA cells
    int limit_slopes; /**< DG2 slope limiter enabled when limit_slopes = ON */

	// adaptation. to merge?
//	NUMERIC_TYPE       xmin;
//	NUMERIC_TYPE       xmax;
//	NUMERIC_TYPE       ymin;
//	NUMERIC_TYPE       ymax;
//	Coordinate xsz;
//	Coordinate ysz;
//	NUMERIC_TYPE       g;
//	NUMERIC_TYPE       time;
//	NUMERIC_TYPE       manning;

};

// Solver settings
struct Solver{
	NUMERIC_TYPE t;
	NUMERIC_TYPE g;
	NUMERIC_TYPE divg;
	NUMERIC_TYPE cfl;
	int    ts_multiple; // channel timestep multiple for running 1D decoupled from 2D
	long   Nit, itCount;
	NUMERIC_TYPE Sim_Time;
	NUMERIC_TYPE InitTstep; // Maximum timestep
	NUMERIC_TYPE Tstep;  // Adapting timestep
	NUMERIC_TYPE MinTstep;  // Stores minimum timestep during simulation
	NUMERIC_TYPE SolverAccuracy;
	int dynsw; // Switch for full dynamic steady state or diffusive steady state
	NUMERIC_TYPE Impfactor;
	NUMERIC_TYPE Hds;
	NUMERIC_TYPE vol1, vol2;
	NUMERIC_TYPE Qerror;
	NUMERIC_TYPE Verror;
	NUMERIC_TYPE FArea; // Store flooded area
	NUMERIC_TYPE DepthThresh, MomentumThresh, MaxHflow;
	NUMERIC_TYPE dhlin;
	NUMERIC_TYPE htol;
	NUMERIC_TYPE Qlimfact; // MT added to allow user to relax Qlimit
	NUMERIC_TYPE itrn_time;
	NUMERIC_TYPE itrn_time_now;
	NUMERIC_TYPE SGCtmpTstep; // JCN added to enable time step calculatin in UpdateH for SGC method
	time_t time_start;
	time_t time_finish;
	time_t time_check;
	NUMERIC_TYPE theta; //GAMA added for q-centred scheme
	int fricSolver2D; //GAMA: Solves the friction term using the vectorial (2D) scheme
	NUMERIC_TYPE maxH; /**< maximum H in the domain at the current time */
	NUMERIC_TYPE krivodonova_threshold; /**< DG2 slope detector */
	NUMERIC_TYPE SpeedThresh; /**< FV1/DG2 threshold for friction application */
    NUMERIC_TYPE DG2DepthThresh; /**< Threshold above which DG2 L1 operator is activated */
    NUMERIC_TYPE DG2ThinDepthTstep; /**< Tstep assigned to cells with thin depths */

	// adaptation
	NUMERIC_TYPE epsilon; // error threshold for adaptation
	int L; // max resolution level for adaptation
};


//-------------------------------------------
// ChannelSegmentType
struct ChannelSegmentType{
	NUMERIC_TYPE *Chandx;
	NUMERIC_TYPE *Shalf;
	NUMERIC_TYPE *Chainage;
	NUMERIC_TYPE *ChanQ; // only for recording Q for output in profile
	NUMERIC_TYPE *A;
	NUMERIC_TYPE *NewA;
	NUMERIC_TYPE *ChanWidth;
	NUMERIC_TYPE *ChanN;
	int   *ChanX;
	int 	*ChanY;
	// BC_Val used in case of fixed flow (otherwise set to 0)
	NUMERIC_TYPE *Q_Val;
	// time series indexed by bci (boundary condition index) //TFD
	// Q_TimeSeries used in case of var e.e. HVAR or QVAR (otherwise set to NULL)
	TimeSeries **Q_TimeSeries;
	ESourceType *Q_Ident;
	char  *Q_Name;
	NUMERIC_TYPE *BankZ;
	int chsz;
	int Next_Segment;
	int Next_Segment_Loc;
	int N_Channel_Segments;
	NUMERIC_TYPE JunctionH; // allows recording of H data for dummy junction point of tributary, without overwriting main channel info
	NUMERIC_TYPE JunctionDEM; // allows recording of DEM data for dummy junction point of tributary, without overwriting main channel info
};

// DamType //FEOL
// structure is needed for DynamicEdges as list per dam is needed
struct DamEdge {
	int *EdgeCell;
};
struct DamData {
	NUMERIC_TYPE *DynamicEdgeData;
	int *DynamicEdge;
	//NUMERIC_TYPE *Outputcell; Not needed
	NUMERIC_TYPE *Volmax;
	NUMERIC_TYPE *DamArea;
	NUMERIC_TYPE *InitialHeight; //Initial Height of Dam
	NUMERIC_TYPE *DamHeight; // Internal Depth of Dam for Vol calc
	NUMERIC_TYPE *SpillWidth;
	NUMERIC_TYPE *Spill_Cd;
	NUMERIC_TYPE *SpillHeight; // Crest Height of Spill
	NUMERIC_TYPE *DamOperationQ;
	int *DamOperationCode; //1 is default and removes MeanQ for Dam as Operation Rule
	NUMERIC_TYPE *DamMeanQ;
	NUMERIC_TYPE *DamVin;
	NUMERIC_TYPE *SpillQ;
	NUMERIC_TYPE *DamTotalQ;
	NUMERIC_TYPE *DamVol;
	NUMERIC_TYPE DamLoss;
	NUMERIC_TYPE *AnnualRelease;
	NUMERIC_TYPE *OP7_Kappa;
	NUMERIC_TYPE *OutputCellX;
	NUMERIC_TYPE *OutputCellY;
	int NumDams;
	int *Edgenos;
	int TotalEdge;
	int DamMaxH;//Check FEOL 21July
	int *DamYear;
	};

//-------------------------------------------
/* QID7_Store // CCS for temp storage of trib boundary condition info when (Q_Ident_tmp[i]==7) in LoadRiver. A vector containing these
structures is built in LoadRiver function and used in UpdateChannelsVector function. */
struct QID7_Store{
	int trib;
	int Next_Segment_Loc;
	int chseg;
	int RiverID;
};
//-------------------------------------------


/*
*****************************************************************************

Define the function prototypes
---------------------
Prototypes split into approximate groups that correspond to file locations
for easier editing.

*****************************************************************************
*/

// Input prototypes - input.cpp
void ReadConfiguration
(
	int argc,
	char *argv[],
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	const int verbose
);

int ReadVerboseMode(int argc, char *argv[]);
void ReadCommandLine(int argc, char *argv[], Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose);
void ReadParamFile(char *, Fnames *, States *, Pars *, Solver*, int);
void CheckParams(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, int verbose);
void LoadDEM(Fnames *, States *, Pars *, Arrays *, const int verbose);
FILE* LoadDomainGeometry(const char* filename, Pars *Parptr, const int verbose, NUMERIC_TYPE& no_data_value);
void LoadDEMData(Pars*, NUMERIC_TYPE *DEM, FILE *fp, NUMERIC_TYPE file_nodata_value);
void LoadManningsn(Fnames *, Pars *, Arrays *, const int verbose);
void LoadDistInfil(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, const int verbose);
void LoadSGCManningsn(Fnames *, Pars *, Arrays *, const int verbose);
void LoadSGCdirn(Fnames *, Pars *, Arrays *, const int verbose);
void LoadRiverNetwork(Fnames *, States *, Pars *, vector<ChannelSegmentType> *, Arrays *, vector<QID7_Store> *, vector<int> *, const int verbose); // CCS
void LoadRiver(Fnames *, States *, Pars *, vector<ChannelSegmentType> *, Arrays *, vector<QID7_Store> *, vector<int> *, const int verbose);
void UpdateChannelsVector(States *, ChannelSegmentType *, vector<QID7_Store> *, QID7_Store *, int *); // CCS
void LoadStart(Fnames *, States *, Pars *, Arrays *, SGCprams *, const int verbose);
void LoadStartQ2D(Fnames*, Pars*, Arrays*, const int verbose);
void LoadBCs(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, BoundCs *BCptr, const int verbose);
void LoadBCVar(Fnames *, States *, Pars *, BoundCs *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *, const int verbose);
void LoadWeir(Fnames *, States *, Pars *, Arrays *, const int verbose);
void LoadStages(Fnames *, States *, Pars *, Stage *, const int verbose);
void LoadGauges(Fnames *, States *, Pars *, Stage *, const int verbose);
void LoadPor(Fnames *, States *, Pars *, Arrays *, const int verbose);
void LoadEvap(Fnames *, Arrays *, const int verbose);
void LoadRain(Fnames *, Arrays *, const int verbose);
void LoadRainmask(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, States *Statesptr, const int verbose);
void LoadSGC(Fnames *Fnameptr, Pars *Parptr, Arrays *Arrptr, States *Statesptr, const int verbose);
void LoadBinaryStart(Fnames *, States *, Pars *, Arrays *, SGCprams *SGCptr, const int verbose);
void LoadSGCChanPrams(Fnames *, States *, Pars *, SGCprams *, const int verbose);
void LoadDamPrams(Fnames *Fnameptr, States *Statesptr, Pars *Parptr, DamData *Damptr, const int verbose); //FEOL
void LoadDamMask(Fnames *Fnameptr,Pars *Parptr, Arrays *Arrptr, DamData *Damptr, const int verbose);  //FEOL
FILE* fopen_or_die(const char * filename, const char* mode, const char* message = "", const int verbose = OFF);

// LISFLOOD Solution prototypes - iterateq.cpp
void IterateQ(Fnames *, Files *, States *, Pars *, Solver*, BoundCs *, Stage *, ChannelSegmentType *, Arrays *, SGCprams *, vector<int> *, int *, vector<ChannelSegmentType> *, const int verbose);
void UpdateH(States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *);

// Floodplain prototypes - fp_flow.cpp
void FloodplainQ(States *, Pars *, Solver *, Arrays *, SGCprams *);
NUMERIC_TYPE CalcFPQx(int i, int j, States *, Pars *, Solver *, Arrays *, NUMERIC_TYPE * TSptr);
NUMERIC_TYPE CalcFPQy(int i, int j, States *, Pars *, Solver *, Arrays *, NUMERIC_TYPE * TSptr);
int MaskTest(int m1, int m2);
int MaskTestAcc(int m1);

// Channel prototypes - ch_flow.cpp
void SetChannelStartH(States *Statesptr, Pars *Parptr, Arrays *Arrptr, ChannelSegmentType *ChannelSegments, vector<int> *, int *);
void SetChannelStartHfromQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr, ChannelSegmentType *ChannelSegments, Solver *, vector<int> *, int *);
void CalcChannelStartQ(States *Statesptr, Pars *Parptr, Arrays *Arrptr, ChannelSegmentType *ChannelSegments, vector<int> *, int *);
void ChannelQ(NUMERIC_TYPE deltaT, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, vector<int> *, int *);
NUMERIC_TYPE CalcA(NUMERIC_TYPE n, NUMERIC_TYPE s, NUMERIC_TYPE w, NUMERIC_TYPE Q);
NUMERIC_TYPE BankQ(int chani, ChannelSegmentType *, Pars *, Arrays *);
NUMERIC_TYPE ChannelVol(States *, Pars *, ChannelSegmentType *, Arrays *);
NUMERIC_TYPE CalcQ(NUMERIC_TYPE n, NUMERIC_TYPE s, NUMERIC_TYPE w, NUMERIC_TYPE h);
NUMERIC_TYPE Newton_Raphson(NUMERIC_TYPE Ai, NUMERIC_TYPE dx, NUMERIC_TYPE a0, NUMERIC_TYPE a1, NUMERIC_TYPE c, Solver *);

// Diffusive channel solver specific functions
void ChannelQ_Diff(NUMERIC_TYPE deltaT, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, vector<int> *, int *);
void bandec(NUMERIC_TYPE **a, int n, int m1, int m2, NUMERIC_TYPE **al, int indx[], NUMERIC_TYPE &d);
void banbks(NUMERIC_TYPE **a, int n, int m1, int m2, NUMERIC_TYPE **al, int indx[], NUMERIC_TYPE b[]);
void SWAP(NUMERIC_TYPE &a, NUMERIC_TYPE &b);
void calcF(NUMERIC_TYPE *x, NUMERIC_TYPE *xn, NUMERIC_TYPE *f, NUMERIC_TYPE dt, ChannelSegmentType *csp, Pars *Parptr, Arrays *Arrptr, NUMERIC_TYPE Qin, int chseg, NUMERIC_TYPE WSout, int HoutFREE, Solver *Solverptr, int low);
void calcJ(NUMERIC_TYPE *x, NUMERIC_TYPE *xn, NUMERIC_TYPE **J, NUMERIC_TYPE dt, ChannelSegmentType *csp, int chseg, int HoutFREE);
NUMERIC_TYPE norm(NUMERIC_TYPE *x, int n);
NUMERIC_TYPE CalcEnergySlope(NUMERIC_TYPE n, NUMERIC_TYPE w, NUMERIC_TYPE h, NUMERIC_TYPE Q);
//void precond(NUMERIC_TYPE **a, int n);

// Acceleration floodplain solver
NUMERIC_TYPE CalcFPQxAcc(int i, int j, States *, Pars *, Solver *, Arrays *);
NUMERIC_TYPE CalcFPQyAcc(int i, int j, States *, Pars *, Solver *, Arrays *);
NUMERIC_TYPE CalcMaxH(Pars *, Arrays *);
void CalcT(Pars *, Solver *, Arrays *);
void UpdateQs(Pars *, Arrays *);




// Boundary prototypes - boundary.cpp
void BCs(States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *);
void BoundaryFlux(States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *);
NUMERIC_TYPE InterpolateTimeSeries(TimeSeries *timeSeries, NUMERIC_TYPE t);

/*
 * Implements a free boundary condition for irregular domains by zeroing
 * water depths over cells with DEM value nodata_elevation
 */
void drain_nodata_water(Pars*, Solver*, BoundCs*, Arrays*);

NUMERIC_TYPE RoeBCy(int edge, int p0, int p1, int pq0, NUMERIC_TYPE z0, NUMERIC_TYPE z1, NUMERIC_TYPE hl, NUMERIC_TYPE hr, NUMERIC_TYPE hul, NUMERIC_TYPE hur, NUMERIC_TYPE hvl, NUMERIC_TYPE hvr, States *Statesptr, Pars *Parptr, Solver *Solverptr, Arrays *Arrptr);
NUMERIC_TYPE RoeBCx(int edge, int p0, int p1, int pq0, NUMERIC_TYPE z0, NUMERIC_TYPE z1, NUMERIC_TYPE hl, NUMERIC_TYPE hr, NUMERIC_TYPE hul, NUMERIC_TYPE hur, NUMERIC_TYPE hvl, NUMERIC_TYPE hvr, States *Statesptr, Pars *Parptr, Solver *Solverptr, Arrays *Arrptr);


// Optional addon protoypes
// chkpnt.cpp
void ReadCheckpoint(Fnames *, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, const int verbose);
void WriteCheckpoint(Fnames *, States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *, const int verbose);
// infevap.cpp
void FPInfiltration(Pars *, Solver *, Arrays *);
void Evaporation(Pars *, Solver *, Arrays *);
void Rainfall(Pars *, Solver *, Arrays *);
void FlowDirDEM(Pars *, Arrays *, States *, BoundCs *); // Calculate routing intervals and flow directions from DEM for rainfall component CCS 13/03/2012
void Routing(States *, Pars *, Solver *, Arrays *); // Route shallow flows from rainfall CCS 14/03/2012
// por_flow.cpp
NUMERIC_TYPE CalcFPQxPor(int i, int j, States *, Pars *, Solver *, Arrays *);
NUMERIC_TYPE CalcFPQyPor(int i, int j, States *, Pars *, Solver *, Arrays *);
NUMERIC_TYPE PorArea(int t, int j, Pars *, Arrays *);
// weir_flow.cpp
NUMERIC_TYPE CalcWeirQx(int i, int j, const Pars *, const Arrays *, const Solver *, const States *, const SGCprams *);
NUMERIC_TYPE CalcWeirQy(int i, int j, const Pars *, const Arrays *, const Solver *, const States *, const SGCprams *);

// Utility prototypes - util.cpp
NUMERIC_TYPE DomainVol(States *, Pars *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *);
void SmoothBanks(Pars *, Solver *, ChannelSegmentType *, Arrays *, vector<ChannelSegmentType> *, const int verbose);
NUMERIC_TYPE x_centre(Pars *Parptr, const int i);
NUMERIC_TYPE y_centre(Pars *Parptr, const int j);
NUMERIC_TYPE x_vertex(Pars *Parptr, const int i);
NUMERIC_TYPE y_vertex(Pars *Parptr, const int j);
void DryCheck(Pars *, Solver *, Arrays *);
int signR(NUMERIC_TYPE a);
void UpdateV(States *, Pars *, Solver *, BoundCs *, ChannelSegmentType *, Arrays *);
void InitFloodplainQ(States *, Pars *, Solver *, Arrays *);
//NUMERIC_TYPE CalcVirtualGauge(int i, Pars *, Arrays *, Stage *);
NUMERIC_TYPE CalcVirtualGauge(const int gauge_i, const int grid_cols_padded,
	const NUMERIC_TYPE * qx_grid, const NUMERIC_TYPE * qy_grid,
	Stage *Locptr);

void CalcArrayDims(States *, Pars *, Arrays *);

// Ouput prototypes - output.cpp
void fileoutput(Fnames *, States *, Pars *, Arrays *);
void write_regular_output(Fnames *, Solver *, States *, Pars *, Arrays *, SGCprams *);

// MT new general purpose functions
void write_ascfile(const char *root, int SaveNumber, const char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *dem, int outflag, States *, Pars *); // general purpose ascii write routine
void write_ascfile(const char *root, int SaveNumber, const char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *dem, int outflag, States *Statesptr, Pars *Parptr, NUMERIC_TYPE depth_thresh);
void write_ascfile(const char *root, int SaveNumber, const char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *dem, int outflag, States *Statesptr, Pars *Parptr, NUMERIC_TYPE depth_thresh, const char* format_specifier);
void write_ascfileDScaled(const char* root, int SaveNumber, const char* extension, NUMERIC_TYPE* data, NUMERIC_TYPE* dem, int outflag, States* Statesptr, Pars* Parptr, NUMERIC_TYPE depth_thresh, const char* format_specifier);
void write_binrasterfile(const char *root, int SaveNumber, const char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *dem, int outflag, States *, Pars *); // general purpose binary write routine
void write_binrasterfile(const char *root, int SaveNumber, const char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *dem, int outflag, States *, Pars *, NUMERIC_TYPE depth_thresh); // general purpose binary write routine
void write_ascfile_SGCf(char *root, int SaveNumber, char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *SGC_BankFullHeight_grid, States *, Pars *); // specific routine for SGC floodplain depth export ascii
void write_binrasterfile_SGCf(char *root, int SaveNumber, char *extension, NUMERIC_TYPE *data, NUMERIC_TYPE *SGC_BankFullHeight_grid, States *, Pars *); // specific routine for SGC floodplain depth export binary
void write_profile(char *root, int SaveNumber, char *extension, States *Statesptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, Pars *Parptr, vector<int> *, int *RiversIndexPtr); // write river channel profiles 
void debugfileoutput(Fnames *, States *, Pars *, Arrays *); // Debug option file output (currently the modified DEM and the channel and trib masks
void printversion(int verbose); // output program version header
int fexist(char *filename); // check if file exists

// TRENT functions
NUMERIC_TYPE maximum(NUMERIC_TYPE a, NUMERIC_TYPE b, NUMERIC_TYPE c);
NUMERIC_TYPE CalcFPQxRoe(int i, int j, States *, Pars *, Solver *, Arrays *);
NUMERIC_TYPE CalcFPQyRoe(int i, int j, States *, Pars *, Solver *, Arrays *);
void UpdateQsRoe(Pars *, Solver *, Arrays *);
//void UpdateHRoe(States *, Pars *, Solver *,BoundCs *,ChannelSegmentType *,Arrays *);
void CalcTRoe(Pars *, Solver *, Arrays *);
