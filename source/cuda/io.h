#pragma once
#include <future>
#include "../geometry.h"
#include "../lisflood.h"

namespace lis
{

struct AsciiRaster
{
	static void read_header
	(
		FILE* file,
		Geometry& geometry,
		NUMERIC_TYPE& no_data_value
	);

	static void read
	(
		FILE* file,
		NUMERIC_TYPE* array,
		Geometry& geometry,
		int pitch,
		int offset
	);

	static void replace_no_data
	(
		NUMERIC_TYPE* array,
		Geometry& geometry,
		int pitch,
		int offset,
		NUMERIC_TYPE original_no_data_value,
		NUMERIC_TYPE new_value
	);

	template<typename F>
	static void write
	(
		FILE* file,
		F array,
		Geometry& geometry,
		int pitch,
		int offset,
		NUMERIC_TYPE no_data_value,
		int outflag,
		int precision = DEFAULT_PRECISION
		
	);

	static void match_cell_dimensions_or_die
	(
		Geometry& expected,
		Geometry& candidate,
		const char* message
	);
};

struct StartFile
{
	static void load
	(
		const char* filename,
		NUMERIC_TYPE* array,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	static void subtract_dem
	(	
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* DEM,
		Geometry& geometry,
		int pitch,
		int offset
	);

private:
	static void zero_no_data_values
	(
		NUMERIC_TYPE* array,
		Geometry& geometry,
		int pitch,
		int offset,
		NUMERIC_TYPE no_data_value
	);
};

struct Snapshot
{
	static std::future<void> initialise_async();

	template<typename F>
	static std::future<void> write_async
	(
		F array,
		Geometry& geometry,
		int pitch,
		int offset,
		const char* prefix,
		int counter,
		const char* suffix,
		int verbose,
		int outflag,
		const int& call_gzip,
		int precision = DEFAULT_PRECISION,
		NUMERIC_TYPE no_data_value = NULLVAL
	);

	template<typename F>
	static void write
	(
		F array,
		Geometry& geometry,
		int pitch,
		int offset,
		const char* prefix,
		int counter,
		const char* suffix,
		int verbose,
		int outflag,
		const int& call_gzip,
		int precision = DEFAULT_PRECISION,
		NUMERIC_TYPE no_data_value = NULLVAL
		
	);

	template<typename F>
	static void write_max
	(
		F array,
		Geometry& geometry,
		int pitch,
		int offset,
		const char* prefix,
		const char* suffix,
		int verbose,
		int outflag,
		int precision = DEFAULT_PRECISION,
		NUMERIC_TYPE no_data_value = NULLVAL

	);
};

}

#include "io.templates.cpp"
