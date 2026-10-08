/*
* file_tool.h
*
*  Created on: 16 Jun 2014
*      Author: td14281
* Toby Dunne - read/write file in ascfile format based on input (no dependency on globals)
*/
#pragma once

void write_time(const int SaveNumber, const NUMERIC_TYPE curr_time,
	const OutputParams* output_params);

#define NETCDF_IGNORE -1

#if _NUMERIC_MODE == 1
#define NC_GET_VARA_NUMERIC_TYPE nc_get_vara_double
#else
#define NC_GET_VARA_NUMERIC_TYPE nc_get_vara_float
#endif

void StartNetCDF(int grid_cols, int grid_rows,
	char* filename,
	NetCDFState * netcdf_state,
	const int save_depth,
	const int save_elev,
	const int save_Qs,
	const int save_Velocity,
	const int save_sub_grid_Velocity,
	const int save_sub_grid);

void CloseNetCDF(NetCDFState * netcdf_state);

void CloseNetCDF(int ncid);

/*
writes grid - data not padded
netcdf_varid: a var_id defined in netcdf_state,  NETCDF_IGNORE if the grid is not written to netcdf

*/
void write_grid(const char *root, const int SaveNumber,
	int netcdf_varid,
	const char *extension,
	const NUMERIC_TYPE *data,
	const int grid_cols, const int grid_rows, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params);

/*
writes grid - data may or may not be padded.
            - if non-padded data required, padding will be removed using tmp_grid

netcdf_varid: a var_id defined in netcdf_state,  NETCDF_IGNORE if the grid is not written to netcdf
*/
void write_grid(const char *root, const int SaveNumber,
	int netcdf_varid,
	const char *extension,
	const NUMERIC_TYPE *data, NUMERIC_TYPE *tmp_grid,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params);

/*
netcdf_varid: a var_id defined in netcdf_state,  NETCDF_IGNORE if the grid is not written to netcdf
*/
void write_sum_grid(const char *root, const int SaveNumber, 
	int netcdf_varid,
	const char *extension, const NUMERIC_TYPE * dataA, const NUMERIC_TYPE * dataB, NUMERIC_TYPE * tmp_grid, const NUMERIC_TYPE a_depth_thresh,
	const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size, OutputParams* output_params);

void write_file_asc(const char *filename, const NUMERIC_TYPE *data, const int grid_cols, const int grid_rows, const int grid_cols_padded, const NUMERIC_TYPE xllcorner, const NUMERIC_TYPE yllcorner, const NUMERIC_TYPE cell_size);

void read_file_asc(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size);

void read_file_bin(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size);

/// 
/// file type is detected from contents (ascii or binary)
///
void read_file(const char * filename, const NUMERIC_TYPE nodata, int * num_cols, int * num_rows, NUMERIC_TYPE** data, NUMERIC_TYPE * xllcorner, NUMERIC_TYPE * yllcorner, NUMERIC_TYPE * cell_size);

void compare_grids(char* dir1, char* dir2, char* prefix, char* ext);

bool read_file_netCDF(NetCDFVariable *ncvar, NUMERIC_TYPE curr_time);
void read_file_netCDF_start(const char *ncfilename, const char *varname, NetCDFVariable *ncvar);
