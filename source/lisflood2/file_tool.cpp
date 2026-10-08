#include "../lisflood.h"
#include "../VersionHistory.h"
#include "lisflood2.h"
#include "file_tool.h"
#include "../utility.h"

#if _NETCDF == 1
#include <netcdf.h>
#endif

/* Handle errors by printing an error message and exiting with a
* non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#if _NUMERIC_MODE == 1
#define nc_put_att_NUMERIC_TYPE nc_put_att_double
#define nc_put_var1_NUMERIC_TYPE nc_put_var1_double
#define NC_NUMERIC_TYPE NC_DOUBLE

#else
#define nc_put_att_NUMERIC_TYPE nc_put_att_float
#define nc_put_var1_NUMERIC_TYPE nc_put_var1_float
#define NC_NUMERIC_TYPE NC_FLOAT
#endif

#define UNITS_LABEL "units"
#define LONG_NAME_LABEL "long_name"
#define MISSING_VALUE_LABEL	"missing_value"

#define HEIGHT_UNITS_VALUE "meters" 
#define FLOW_UNITS_VALUE "meters3/second" 
#define VELOCITY_UNITS_VALUE "meters/second" 
#define TIME_DURATION_HOURS_UNITS_VALUE "hours"

#if _NETCDF == 1
/*
  init grid with 2 spacial dimentions and one time dimention
  */
void NetCDF_init_3d_grid(int dimid_y, int dimid_x, char* name, char* long_name, char* units,
	NetCDFState * netcdf_state,
	int * out_var_id)
{
	/* Set chunking, shuffle, and deflate. */
	int shuffle = NC_SHUFFLE;
	int deflate = 1;
	int deflate_level = 1;
	NUMERIC_TYPE no_data_value = NULLVAL;
	int retval;
	int dimids[3];
	dimids[0] = netcdf_state->dimid_time; // unlimited dimention must be first
	dimids[1] = dimid_y; // y dimention first due to how grid is stored
	dimids[2] = dimid_x;

	if ((retval = nc_def_var(netcdf_state->ncid, name, NC_NUMERIC_TYPE, 3,
		dimids, out_var_id)))
		ERR(retval);

	int varid = *out_var_id;

	if (long_name != NULL)
	{
		if ((retval = nc_put_att_text(netcdf_state->ncid, varid, LONG_NAME_LABEL,
			strlen(long_name), long_name)))
			ERR(retval);
	}

	if ((retval = nc_put_att_text(netcdf_state->ncid, varid, UNITS_LABEL,
		strlen(HEIGHT_UNITS_VALUE), HEIGHT_UNITS_VALUE)))
		ERR(retval);

	if (retval = nc_put_att_NUMERIC_TYPE(netcdf_state->ncid, varid, MISSING_VALUE_LABEL,
		NC_NUMERIC_TYPE, 1, &no_data_value))
		ERR(retval);

	// don't specify chunking - default used
	//if ((retval = nc_def_var_chunking(netcdf_state->ncid, varid, 0, &chunks[0])))
	//	ERR(retval);

	if ((retval = nc_def_var_deflate(netcdf_state->ncid, varid, shuffle, deflate,
		deflate_level)))
		ERR(retval);

}

/*
init grid with 2 spacial dimentions
*/
void NetCDF_init_2d_grid(int dimid_y, int dimid_x, char* name, char* longname, char* units,
	NetCDFState * netcdf_state,
	int * out_var_id)
{
	/* Set chunking, shuffle, and deflate. */
	int shuffle = NC_SHUFFLE;
	int deflate = 1;
	int deflate_level = 1;
	NUMERIC_TYPE no_data_value = NULLVAL;
	int retval;
	int dimids[2];
	dimids[0] = dimid_y; // y dimention first due to how grid is stored
	dimids[1] = dimid_x;

	if ((retval = nc_def_var(netcdf_state->ncid, name, NC_NUMERIC_TYPE, 2,
		dimids, out_var_id)))
		ERR(retval);

	int varid = *out_var_id;

	if (longname != NULL)
	{
		if ((retval = nc_put_att_text(netcdf_state->ncid, varid, LONG_NAME_LABEL,
			strlen(longname), longname)))
			ERR(retval);
	}

	if ((retval = nc_put_att_text(netcdf_state->ncid, varid, UNITS_LABEL,
		strlen(HEIGHT_UNITS_VALUE), HEIGHT_UNITS_VALUE)))
		ERR(retval);

	if (retval = nc_put_att_NUMERIC_TYPE(netcdf_state->ncid, varid, MISSING_VALUE_LABEL,
		NC_NUMERIC_TYPE, 1, &no_data_value))
		ERR(retval);

	// don't specify chunking - default used
	//if ((retval = nc_def_var_chunking(netcdf_state->ncid, varid, 0, &chunks[0])))
	//	ERR(retval);

	if ((retval = nc_def_var_deflate(netcdf_state->ncid, varid, shuffle, deflate,
		deflate_level)))
		ERR(retval);

}

void StartNetCDF_impl(int grid_cols, int grid_rows,
	char* filename,
	NetCDFState * netcdf_state,
	const int save_depth,
	const int save_elev,
	const int save_Qs,
	const int save_Velocity,
	const int save_sub_grid_Velocity,
	const int save_sub_grid)
{
	netcdf_state->init_done = ON;

#define TIME_UNITS_VALUE "seconds since 1970-01-01 00:00:00 0:00"
#define TIME_HOURS_UNITS_VALUE "hours since 1970-01-01 00:00:00 0:00"

	int retval;
	/* Create the file. The NC_NETCDF4 parameter tells netCDF to create
	* a file in netCDF-4/HDF5 standard. */
	if ((retval = nc_create(filename, NC_NETCDF4, &netcdf_state->ncid)))
		ERR(retval);

	char source_string[256];
	sprintf(source_string, "LISFLOOD - FP version %d.%d.%d(%s)\n", LF_VersionMajor, LF_VersionMinor, LF_VersionInc, NUMERIC_TYPE_NAME);
	if ((retval = nc_put_att_text(netcdf_state->ncid, NC_GLOBAL, "Source", strlen(source_string), source_string)))
		ERR(retval);

	//#define CONVENTIONS_VALUE "COARDS, GDV"
	//if ((retval = nc_put_att_text(netcdf_state->ncid, NC_GLOBAL, "Conventions", strlen(CONVENTIONS_VALUE), CONVENTIONS_VALUE)))
	//	ERR(retval);

	if ((retval = nc_def_dim(netcdf_state->ncid, "time", NC_UNLIMITED, &netcdf_state->dimid_time)))
		ERR(retval);

	int dimids[1];
	size_t chunking[1];
	dimids[0] = netcdf_state->dimid_time;

	/* Define the time variable. */
	if ((retval = nc_def_var(netcdf_state->ncid, "time", NC_NUMERIC_TYPE, 1,
		dimids, &netcdf_state->varid_time)))
		ERR(retval);

	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_time, UNITS_LABEL,
		strlen(TIME_UNITS_VALUE), TIME_UNITS_VALUE)))
		ERR(retval);

	//#define TIME_TIME_STEP_LABEL "time_step" 
	//#define TIME_TIME_STEP_VALUE "seconds" 
	//	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_time, TIME_TIME_STEP_LABEL,
	//		strlen(TIME_TIME_STEP_VALUE), TIME_TIME_STEP_VALUE)))
	//		ERR(retval);

	chunking[0] = 1;

	if ((retval = nc_def_var_chunking(netcdf_state->ncid, netcdf_state->varid_time, 0, chunking)))
		ERR(retval);

#ifdef LATLONGMODE
	if ((retval = nc_def_dim(netcdf_state->ncid, "lat", grid_cols, &netcdf_state->dimid_x)))
		ERR(retval);
	if ((retval = nc_def_dim(netcdf_state->ncid, "lon", grid_rows, &netcdf_state->dimid_y)))
		ERR(retval);


#define LONG_NAME_Longitude "Longitude"
#define UNITS_Longitude "degrees_east"

#define LONG_NAME_Latitude "Latitude"
#define UNITS_Latitude "degrees_north"

	if ((retval = nc_def_var(netcdf_state->ncid, "lat", NC_NUMERIC_TYPE, 1,
		dimids, &netcdf_state->varid_y)))
		ERR(retval);

	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_y, LONG_NAME_LABEL,
		strlen(LONG_NAME_Latitude), LONG_NAME_Latitude)))
		ERR(retval);
	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_y, UNITS_LABEL,
		strlen(UNITS_Latitude), UNITS_Latitude)))
		ERR(retval);


	if ((retval = nc_def_var(netcdf_state->ncid, "lon", NC_NUMERIC_TYPE, 1,
		dimids, &netcdf_state->varid_x)))
		ERR(retval);

	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_x, LONG_NAME_LABEL,
		strlen(LONG_NAME_Longitude), LONG_NAME_Longitude)))
		ERR(retval);
	if ((retval = nc_put_att_text(netcdf_state->ncid, netcdf_state->varid_x, UNITS_LABEL,
		strlen(UNITS_Longitude), UNITS_Longitude)))
		ERR(retval);



#else

	if ((retval = nc_def_dim(netcdf_state->ncid, "y", grid_rows, &netcdf_state->dimid_y)))
		ERR(retval);

	if ((retval = nc_def_dim(netcdf_state->ncid, "x", grid_cols, &netcdf_state->dimid_x)))
		ERR(retval);

	if ((retval = nc_def_dim(netcdf_state->ncid, "y_edge", grid_rows + 1, &netcdf_state->dimid_y_edge)))
		ERR(retval);

	if ((retval = nc_def_dim(netcdf_state->ncid, "x_edge", grid_cols + 1, &netcdf_state->dimid_x_edge)))
		ERR(retval);

#endif

	if (save_depth == ON)
	{
		NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "depth", NULL, HEIGHT_UNITS_VALUE, netcdf_state, &netcdf_state->varid_depth);
	}

	if (save_elev == ON)
	{
		NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "elevation", NULL, HEIGHT_UNITS_VALUE, netcdf_state, &netcdf_state->varid_depth);
	}

	if (save_Qs == ON)
	{
		NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x_edge, "Qx", "flood plain flow between cells in x direction", FLOW_UNITS_VALUE, netcdf_state, &netcdf_state->varid_qx);
		NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "Qy", "flood plain flow between cells in y direction", FLOW_UNITS_VALUE, netcdf_state, &netcdf_state->varid_qy);

		if (save_sub_grid == ON)
		{
			NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x_edge, "Qcx", "sub grid channel flow between cells in x direction", FLOW_UNITS_VALUE, netcdf_state, &netcdf_state->varid_qcx);
			NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "Qcy", "sub grid channel flow between cells in y direction", FLOW_UNITS_VALUE, netcdf_state, &netcdf_state->varid_qcy);
		}
	}

	if (save_Velocity == ON)
	{
		NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x_edge, "Vx", "flood plain velocity between cells in x direction", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_Vx);
		NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "Vy", "flood plain velocity between cells in y direction", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_Vy);
		//NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "Vc", "flood plain channel cell velocity", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_sgc_Vc);
	}
	if (save_sub_grid == ON && save_sub_grid_Velocity == ON)
	{
		NetCDF_init_3d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x_edge, "SGCVx", "sub grid channel velocity between cells in x direction", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_sgc_Vx);
		NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "SGCVy", "sub grid channel velocity between cells in y direction", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_sgc_Vy);
		NetCDF_init_3d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "SGCVc", "sub grid channel cell velocity", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_sgc_Vc);
	}

	{
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "inittm", "time of initial flood inundation", TIME_HOURS_UNITS_VALUE, netcdf_state, &netcdf_state->varid_inittm);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "totaltm", "total inundation time", TIME_DURATION_HOURS_UNITS_VALUE, netcdf_state, &netcdf_state->varid_totaltm);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "max", "maximum depth", HEIGHT_UNITS_VALUE, netcdf_state, &netcdf_state->varid_max);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "mxe", "maximum elevation", HEIGHT_UNITS_VALUE, netcdf_state, &netcdf_state->varid_mxe);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "maxtm", "time of maximum depth", TIME_HOURS_UNITS_VALUE, netcdf_state, &netcdf_state->varid_maxtm);

		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x_edge, "maxVx", "maximum flood plain velocity between cells in x direction", TIME_HOURS_UNITS_VALUE, netcdf_state, &netcdf_state->varid_maxVx);
		NetCDF_init_2d_grid(netcdf_state->dimid_y_edge, netcdf_state->dimid_x, "maxVy", "maximum flood plain velocity between cells in y direction", TIME_HOURS_UNITS_VALUE, netcdf_state, &netcdf_state->varid_maxVy);

		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "maxVc", "max cell velocity", VELOCITY_UNITS_VALUE, netcdf_state, &netcdf_state->varid_maxVc);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "maxVcd", "depth at max max cell velocity", HEIGHT_UNITS_VALUE, netcdf_state, &netcdf_state->varid_maxVcd);
		NetCDF_init_2d_grid(netcdf_state->dimid_y, netcdf_state->dimid_x, "maxHaz", "maximum hazard", "unknown", netcdf_state, &netcdf_state->varid_maxHaz);
	}


#ifdef LATLONGMODE
	// todo label both axis
#else
	// todo label both axis ??
#endif

}

void CloseNetCDF(NetCDFState * netcdf_state)
{
    CloseNetCDF(netcdf_state->ncid);
}

void CloseNetCDF(int ncid)
{
	int retval;
	if (retval = nc_close(ncid)) ERR(retval);
}

void read_file_netCDF_start(const char *ncfilename, const char *varname,
        NetCDFVariable *ncvar)
{
  int recid, xid, yid;
  int retval;
  
  if ((retval = nc_open(ncfilename, NC_NOWRITE, &ncvar->ncid)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncvar->ncid, varname, &ncvar->varid)))
    ERR(retval);
  if ((retval = nc_inq_varid(ncvar->ncid, "time", &recid)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncvar->ncid, "x", &xid)))
    ERR(retval);
  if ((retval = nc_inq_dimid(ncvar->ncid, "y", &yid)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncvar->ncid, recid, &ncvar->tlen)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncvar->ncid, xid, &ncvar->xlen)))
    ERR(retval);
  if ((retval = nc_inq_dimlen(ncvar->ncid, yid, &ncvar->ylen)))
    ERR(retval);

  size_t start[] = {0};
  size_t tcount[] = {ncvar->tlen};
  ncvar->times = (NUMERIC_TYPE*) calloc(ncvar->tlen, sizeof(NUMERIC_TYPE));
  if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, recid, start, tcount,
              ncvar->times))
    ERR(retval);
  ncvar->time_idx = -1;

  size_t xcount[] = {ncvar->xlen};
  ncvar->xs = (NUMERIC_TYPE*) calloc(ncvar->xlen, sizeof(NUMERIC_TYPE));
  if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, xid, start, xcount,
              ncvar->xs))
    ERR(retval);

  size_t ycount[] = {ncvar->ylen};
  ncvar->ys = (NUMERIC_TYPE*) calloc(ncvar->ylen, sizeof(NUMERIC_TYPE));
  if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, yid, start, ycount,
              ncvar->ys))
    ERR(retval);
}

// assume data are ordered left-to-right, top-to-bottom
// return true if new data needed to be read for the curr_time
bool read_file_netCDF(NetCDFVariable *ncvar, NUMERIC_TYPE curr_time)
{
  int retval;

  size_t old_idx = ncvar->time_idx;
  while (static_cast<int>(ncvar->time_idx) < static_cast<int>(ncvar->tlen)-1 &&
          ncvar->times[ncvar->time_idx+1] <= curr_time)
  {
      ncvar->time_idx++;
  }

  if (old_idx != ncvar->time_idx) {
    if (ncvar->time_idx < ncvar->tlen-1)
    {
        ncvar->dt = ncvar->times[ncvar->time_idx+1]
            - ncvar->times[ncvar->time_idx];
    }
    else
    {
        ncvar->dt = ncvar->times[ncvar->time_idx]
            - ncvar->times[ncvar->time_idx-1];
    }

    size_t start[] = {ncvar->time_idx, 0, 0};
    size_t count[] = {1, ncvar->ylen, ncvar->xlen};
    if (retval = NC_GET_VARA_NUMERIC_TYPE(ncvar->ncid, ncvar->varid, start,
                count, ncvar->data)) ERR(retval);
  }

  return old_idx != ncvar->time_idx;
}

/*
netcdf_varid: must be defined in netcdf_state
*/
void write_file_netCDF(NetCDFState * netcdf_state, int netcdf_varid,
	int time_increment,
	const NUMERIC_TYPE *data,
	const int grid_cols, const int grid_rows, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size)
{
	int retval;
	size_t count[3];
	count[0] = 1; //time
	count[1] = grid_rows; // y
	count[2] = grid_cols; // x

	size_t start[3];
	start[0] = time_increment; //time
	start[1] = 0; // y
	start[2] = 0; // x

	if ((retval = nc_put_vara(netcdf_state->ncid, netcdf_varid, start, count, data)))
		ERR(retval);
}

#endif

void StartNetCDF(int grid_cols, int grid_rows,
	char* filename,
	NetCDFState * netcdf_state,
	const int save_depth,
	const int save_elev,
	const int save_Qs,
	const int save_Velocity,
	const int save_sub_grid_Velocity,
	const int save_sub_grid)
{
#if _NETCDF == 1
	StartNetCDF_impl(grid_cols, grid_rows, filename, netcdf_state, save_depth, save_elev, save_Qs, save_Velocity, save_sub_grid_Velocity, save_sub_grid);
#else
	printf("WARNING: NetCDF selected but not enabled in this build\n");
#endif
}

void write_time(const int SaveNumber, const NUMERIC_TYPE curr_time,
	const OutputParams* output_params)
{
#if _NETCDF == 1
	if (output_params->netcdf_out == ON)
	{
		const NetCDFState * netcdf_state = &output_params->netcdf_state;

		size_t start[1];
		start[0] = SaveNumber;

		int retval;
		if (retval = nc_put_var1_NUMERIC_TYPE(netcdf_state->ncid, netcdf_state->varid_time, start, &curr_time))
			ERR(retval);
	}
#endif
}

///
/// if data is NULL, data is allocated within this function and must be freed by user with memory_free_legacy
///
void read_file_asc(FILE *fp, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size)
{
	char dum[800];
	NUMERIC_TYPE no_data_value = NULLVAL;
	int i, j;
	int ncols, nrows;

	fscanf(fp, "%s %i", dum, &ncols);
	fscanf(fp, "%s %i", dum, &nrows);
	fscanf(fp, "%s %" NUM_FMT"", dum, xllcorner);
	fscanf(fp, "%s %" NUM_FMT"", dum, yllcorner);
	fscanf(fp, "%s %" NUM_FMT"", dum, cell_size);
	fscanf(fp, "%s %" NUM_FMT"", dum, &no_data_value);

	*num_cols = ncols;
	*num_rows = nrows;

	if (*data == NULL)
	{
		*data = memory_allocate_numeric_legacy(ncols * nrows);
	}

	NUMERIC_TYPE *data_grid = *data;

	for (j = 0; j < nrows; j++)
	{
		for (i = 0; i < ncols; i++)
		{
			int index = i + j*ncols;
			NUMERIC_TYPE value;
			fscanf(fp, "%" NUM_FMT"", &value);
			if ((int)value == (int)no_data_value)
				value = nodata;// 1e10;
			data_grid[index] = value;
		}
	}
}

void read_file_asc(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size)
{
	FILE *fp;
	fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		printf("File not found %s\n", filename);
		exit(0);
	}
	read_file_asc(fp, nodata, num_cols, num_rows, data, xllcorner, yllcorner, cell_size);
}

void read_file_bin(FILE *fp, const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size)
{
	NUMERIC_TYPE no_data_value = NULLVAL;
	int i, j;
	int ncols, nrows;

	int count;
	int header_int[2];
	NUMERIC_TYPE header_numeric[4];

	count = fread(header_int, sizeof(int), 2, fp);
	if (count != 2)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}
	count = fread(header_numeric, sizeof(NUMERIC_TYPE), 4, fp);
	if (count != 4)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}

	ncols = header_int[0];
	nrows = header_int[1];

	if (fseek(fp, 0, SEEK_END) != 0)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}

	long file_length = ftell(fp);

	if (file_length != 8 + (4 + ncols * nrows) * sizeof(NUMERIC_TYPE))
	{
		printf("File read error: length mismatch, possible incorrect precision. Expect: %s %s\n", NUMERIC_TYPE_NAME, filename);
		exit(0);
	}


	fseek(fp, 8, SEEK_SET);

	*xllcorner = header_numeric[0];
	*yllcorner = header_numeric[1];
	*cell_size = header_numeric[2];
	no_data_value = header_numeric[3];

	*num_cols = ncols;
	*num_rows = nrows;

	if (*data == NULL)
	{
		*data = memory_allocate_numeric_legacy(ncols * nrows);
	}
	NUMERIC_TYPE *data_grid = *data;

	count = fread(*data, sizeof(NUMERIC_TYPE), ncols * nrows, fp);
	if (count != ncols * nrows)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}
	if (no_data_value != nodata)
	{
		for (j = 0; j < nrows; j++)
		{
			for (i = 0; i < ncols; i++)
			{
				int index = i + j * ncols;
				if (data_grid[index] == no_data_value)
					data_grid[index] = nodata;// 1e10;
			}
		}
	}
}
///
/// if data is NULL, data is allocated within this function and must be freed by user with memory_free_legacy
///
void read_file_bin(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size)
{
	FILE *fp;
	fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		printf("File not found %s\n", filename);
		exit(0);
	}
	read_file_bin(fp, filename, nodata, num_cols, num_rows, data, xllcorner, yllcorner, cell_size);
}

void read_file(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size)
{
	FILE *fp;
	int i;

	fp = fopen(filename, "rb");
	if (fp == NULL)
	{
		printf("File not found %s\n", filename);
		exit(0);
	}

	int count;
	char first[4];
	char firstlower[5];
	count = fread(first, 1, 4, fp);
	if (count != 4)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}
	for (i = 0; i < 4; i++)
	{
		firstlower[i] = tolower(first[i]);
	}
	firstlower[4] = 0;
	if (fseek(fp, 0, SEEK_SET) != 0)
	{
		printf("File read error %s\n", filename);
		exit(0);
	}
	//detect asc or bin
	if (STRCMPi("ncol", firstlower) == 0)
	{
		read_file_asc(fp, nodata, num_cols, num_rows, data, xllcorner, yllcorner, cell_size);
	}
	else
	{
		read_file_bin(fp, filename, nodata, num_cols, num_rows, data, xllcorner, yllcorner, cell_size);
	}
	fclose(fp);
}

void write_file_asc(const char *filename, const NUMERIC_TYPE *data,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size)
{
	int i;
	int j;
	FILE *fp;
	// open file
	fp = fopen(filename, "wb");

	// check file opened ok, if not give warning and exit function
	if (fp == NULL)
	{
		printf("Problems writing to file %s\n", filename);
		return;
	}

	fprintf(fp, "ncols         %i\n", grid_cols);
	fprintf(fp, "nrows         %i\n", grid_rows);
	fprintf(fp, "xllcorner     %.15" NUM_FMT"\n", xllcorner);
	fprintf(fp, "yllcorner     %.15" NUM_FMT"\n", yllcorner);
	fprintf(fp, "cellsize      %.15" NUM_FMT"\n", cell_size);
	fprintf(fp, "NODATA_value  %" NUM_FMT"\n", NULLVAL);

	for (j = 0; j < grid_rows; j++)
	{
		for (i = 0; i < grid_cols; i++)
		{
			fprintf(fp, "%.3" NUM_FMT"\t", data[i + j*(grid_cols_padded)]);
		}
		fprintf(fp, "\n");
	}

	// close file
	fclose(fp);
}

void write_file_bin(const char *filename, const NUMERIC_TYPE *data,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size)
{
	int j;
	FILE *fp;
	// open file
	fp = fopen(filename, "wb");

	// check file opened ok, if not give warning and exit function
	if (fp == NULL)
	{
		printf("Problems writing to file %s\n", filename);
		return;
	}

	const NUMERIC_TYPE nullval = NULLVAL;
	//TODO check if need to be consistent with old version - or can add header for double/float
	// float version won't be compatible with old version
	//fwrite(&sizeofNUMERIC_TYPE, sizeof(int), 1, fp);

	fwrite(&grid_cols, sizeof(int), 1, fp);
	fwrite(&grid_rows, sizeof(int), 1, fp);
	fwrite(&xllcorner, sizeof(NUMERIC_TYPE), 1, fp);
	fwrite(&yllcorner, sizeof(NUMERIC_TYPE), 1, fp);
	fwrite(&cell_size, sizeof(NUMERIC_TYPE), 1, fp);
	fwrite(&nullval, sizeof(NUMERIC_TYPE), 1, fp);

	for (j = 0; j < grid_rows; j++)
	{
		int row_start_index = j * grid_cols_padded;
		fwrite(data + row_start_index, sizeof(NUMERIC_TYPE), grid_cols, fp);
	}

	// close file
	fclose(fp);
}

/**
Save grid with each cell the result of the following:
if dataA > depth_thresh: NULLVAL
otherwise: dataA  + dataB
*/
NUMERIC_TYPE sum_value(NUMERIC_TYPE A, NUMERIC_TYPE B, NUMERIC_TYPE depth_thresh)
{
	NUMERIC_TYPE value;
	if (A <= depth_thresh)
		value = NULLVAL;
	else
		value = A + B;
	return value;
}

void make_filename(const char *root, const int SaveNumber, const char *extension, const OutputParams* output_params, char * filename)
{
	// check if there is a savenumber to add and create filename
	if (SaveNumber >= 0 && SaveNumber <= 9999)
		sprintf(filename, "%s-%.4d%s", root, SaveNumber, extension);
	else if (SaveNumber > 9999)
		sprintf(filename, "%s-%d%s", root, SaveNumber, extension);
	else
		sprintf(filename, "%s%s", root, extension);

	if (output_params->standard_extensions == ON)
	{
		if (output_params->binary_out == ON)
			strcat(filename, ".bin");
		else if (output_params->ascii_out == ON)
			strcat(filename, ".asc");
	}
	else
	{
		if (output_params->binary_out == ON)
			strcat(filename, "b");
	}
}

/*
netcdf_varid: a var_id defined in netcdf_state,  -1 if the grid is not written to netcdf
*/
void write_grid(const char *root, const int SaveNumber,
	int netcdf_varid,
	const char *extension,
	const NUMERIC_TYPE *data, const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params)
{
	char filename[1024];
	char tmp_sys_com[1024 + 11];

	int write_netcdf = 0;
#if _NETCDF == 1
	if (output_params->netcdf_out == ON &&
		netcdf_varid != NETCDF_IGNORE && // some grids are not written into netcdf
		grid_cols == grid_cols_padded) // only write non padded data to netcdf
	{
		write_file_netCDF(&output_params->netcdf_state,
			netcdf_varid, SaveNumber,
			data, grid_cols, grid_rows, xllcorner, yllcorner, cell_size);
		write_netcdf = 1;
	}
#endif
	if (output_params->binary_out == ON)
	{
		make_filename(root, SaveNumber, extension, output_params, filename);
		write_file_bin(filename, data, grid_cols, grid_rows, grid_cols_padded, xllcorner, yllcorner, cell_size);

		// check if we need to zip the file up
		if (output_params->call_gzip == ON)
		{
			sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", filename);
			system(tmp_sys_com);
		}
	}
	if (output_params->ascii_out == ON ||
		(output_params->netcdf_out == ON && write_netcdf == 0)) // if netcdf not implemented for this grid, default to ascii
	{
		make_filename(root, SaveNumber, extension, output_params, filename);
		write_file_asc(filename, data, grid_cols, grid_rows, grid_cols_padded, xllcorner, yllcorner, cell_size);

		// check if we need to zip the file up
		if (output_params->call_gzip == ON)
		{
			sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", filename);
			system(tmp_sys_com);
		}
	}
}

void write_grid(const char *root, const int SaveNumber,
	int netcdf_varid, const char *extension,
	const NUMERIC_TYPE *data,
	const int grid_cols, const int grid_rows, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params)
{
	write_grid(root, SaveNumber,
		netcdf_varid, extension,
		data, grid_cols, grid_rows, grid_cols, xllcorner, yllcorner, cell_size, output_params);
}

/*
netcdf_varid: a var_id defined in netcdf_state,  -1 if the grid is not written to netcdf
*/
void write_grid(const char *root, const int SaveNumber,
	int netcdf_varid,
	const char *extension,
	const NUMERIC_TYPE *data, NUMERIC_TYPE *tmp_grid,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params)
{
	int cols = grid_cols_padded;

	// if writing netcdf - need to remove padding
	if (grid_cols != grid_cols_padded &&
		netcdf_varid != NETCDF_IGNORE && output_params->netcdf_out == ON)
	{
		for (int j = 0; j < grid_rows; j++)
		{
			int row_start_padded = j * grid_cols_padded;
			int row_start = j * grid_cols;
			memcpy(tmp_grid + row_start, data + row_start_padded, sizeof(NUMERIC_TYPE)*grid_cols);
		}
		cols = grid_cols;
		write_grid(root, SaveNumber, netcdf_varid, extension, tmp_grid, grid_cols, grid_rows, cols, xllcorner, yllcorner, cell_size, output_params);
	}
	else
	{
		write_grid(root, SaveNumber, netcdf_varid, extension, data, grid_cols, grid_rows, cols, xllcorner, yllcorner, cell_size, output_params);
	}
}

/**
Save grid with each cell the result of the following:
if dataA > depth_thresh: NULLVAL
otherwise: dataA  + dataB
*/
void write_sum_grid(const char *root,
	int netcdf_varid,
	const int SaveNumber, const char *extension, const NUMERIC_TYPE * dataA, const NUMERIC_TYPE * dataB, NUMERIC_TYPE * tmp_grid, const NUMERIC_TYPE depth_thresh,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params)
{
	for (int j = 0; j < grid_rows; j++)
	{
		for (int i = 0; i < grid_cols; i++)
		{
			int index_padded = i + j * grid_cols_padded;
			int index = i + j * grid_cols;

			tmp_grid[index] = sum_value(dataA[index_padded], dataB[index_padded], depth_thresh);
		}
	}

	write_grid(root, SaveNumber, netcdf_varid, extension, tmp_grid, grid_cols, grid_rows, grid_cols, xllcorner, yllcorner, cell_size, output_params);
}

void compare_grids(char* message, const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE * dataA, const NUMERIC_TYPE * dataB)
{
	int counts[16];
	NUMERIC_TYPE thresholds[16];

	memset(counts, 0, sizeof(int) * 16);

	int non_zeroA = 0;
	int non_zeroB = 0;
	thresholds[0] = C(0.0);
	thresholds[1] = C(1e-1);
	thresholds[2] = C(1e-2);
	thresholds[3] = C(1e-3);
	thresholds[4] = C(1e-4);
	thresholds[5] = C(1e-5);
	thresholds[6] = C(1e-6);
	thresholds[7] = C(1e-7);
	thresholds[8] = C(1e-8);
	thresholds[9] = C(1e-9);
	thresholds[10] = C(1e-10);
	thresholds[11] = C(1e-11);
	thresholds[12] = C(1e-12);
	thresholds[13] = C(1e-13);
	thresholds[14] = C(1e-14);
	thresholds[15] = C(1e-15);

	int i, j;

	for (j = 0; j < grid_rows; j++)
	{
		int row_start_index = j * grid_cols_padded;
		for (i = 0; i < grid_cols; i++)
		{
			int index = row_start_index + i;
			if (dataA[index] != NULLVAL && dataA[index] != C(0.0))
			{
				non_zeroA++;
			}
			if (dataB[index] != NULLVAL && dataB[index] != C(0.0))
			{
				non_zeroB++;
			}
			NUMERIC_TYPE delta = FABS(dataA[index] - dataB[index]);
			for (int t = 0; t < 16; t++)
			{
				if (delta > thresholds[t])
					counts[t]++;
			}
		}
	}
	printf("%s\t", message);
	printf("%d\t", non_zeroA);
	printf("%d\t", non_zeroB);
	for (int t = 0; t < 16; t++)
	{
		printf("%d\t", counts[t]);
	}
	printf("\n");
}

void compare_grids(char* dir1, char* dir2, char* prefix, char* ext)
{
	int count = 0;
	char filename[256];
	char file1[512];
	char file2[512];
	const NUMERIC_TYPE nodata = -9999;
	int num_cols1; int num_rows1; NUMERIC_TYPE  xllcorner1; NUMERIC_TYPE  yllcorner1; NUMERIC_TYPE  cell_size1;
	int num_cols2; int num_rows2; NUMERIC_TYPE  xllcorner2; NUMERIC_TYPE  yllcorner2; NUMERIC_TYPE  cell_size2;

	NUMERIC_TYPE* data1 = NULL;
	NUMERIC_TYPE* data2 = NULL;

	printf("Compare prefix: %s suffix %s\n", prefix, ext);

	printf("Compare dir1: %s\n", dir1);
	printf("Compare dir2: %s\n", dir2);

	sprintf(filename, "%s-%04d%s", prefix, count, ext);

	printf("%s\tA(!=0)\tB(!=0)\t0\t1e-1\t1e-2\t1e-3\t1e-4\t1e-5\t1e-6\t1e-7\t1e-8\t1e-9\t1e-10\t1e-11\t1e-12\t1e-13\t1e-14\t1e-15\t\n", filename);

	do
	{
		sprintf(filename, "%s-%04d%s", prefix, count, ext);
		sprintf(file1, "%s" FILE_SEP"%s", dir1, filename);
		sprintf(file2, "%s" FILE_SEP"%s", dir2, filename);

		if (!fexist(file1) || !fexist(file2))
			break;

		read_file(file1, nodata, &num_cols1, &num_rows1, &data1, &xllcorner1, &yllcorner1, &cell_size1);
		read_file(file2, nodata, &num_cols2, &num_rows2, &data2, &xllcorner2, &yllcorner2, &cell_size2);

		compare_grids(filename, num_cols1, num_rows1, num_cols1, data1, data2);
		count++;
	} while (true);
}
