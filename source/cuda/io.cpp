#include "io.h"
#include <future>

void lis::AsciiRaster::read_header
(
	FILE* file,
	Geometry& geometry,
	NUMERIC_TYPE& no_data_value
)
{
	fscanf(file, "%*s %i", &(geometry.xsz));
	fscanf(file, "%*s %i", &(geometry.ysz));
	fscanf(file, "%*s %" NUM_FMT, &(geometry.blx));
	fscanf(file, "%*s %" NUM_FMT, &(geometry.bly));
	fscanf(file, "%*s %" NUM_FMT, &(geometry.dx));
	fscanf(file, "%*s %" NUM_FMT, &no_data_value);

	geometry.dy = geometry.dx;
	geometry.tly = geometry.bly + geometry.ysz*geometry.dy;
}

void lis::AsciiRaster::read
(
	FILE* file,
	NUMERIC_TYPE* array,
	Geometry& geometry,
	int pitch,
	int offset
)
{
	for (int j=0; j<geometry.ysz; j++)
	{
		for (int i=0; i<geometry.xsz; i++)
		{
			fscanf(file, "%" NUM_FMT, &(array[j*pitch + i + offset]));
		}
	}
}

void lis::AsciiRaster::replace_no_data
(
	NUMERIC_TYPE* array,
	Geometry& geometry,
	int pitch,
	int offset,
	NUMERIC_TYPE original_no_data_value,
	NUMERIC_TYPE new_value
)
{
	for (int j=0; j<geometry.ysz; j++)
	{
		for (int i=0; i<geometry.xsz; i++)
		{
			NUMERIC_TYPE& value = array[j*pitch + i + offset];
			if (FABS(value - original_no_data_value) < C(1e-6))
			{
				value = new_value;
			}
		}
	}
}

void lis::AsciiRaster::match_cell_dimensions_or_die
(
	Geometry& expected,
	Geometry& candidate,
	const char* message
)
{
	if (expected.xsz != candidate.xsz || expected.ysz != candidate.ysz)
	{
		fprintf(stderr, "ERROR: %s: expected %ix%i cells but got %ix%i. "
				"Aborting.", message, expected.xsz, expected.ysz,
				candidate.xsz, candidate.ysz);
		exit(1);
	}
}

void lis::StartFile::load
(
	const char* filename,
	NUMERIC_TYPE* array,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	FILE* file = fopen_or_die(filename, "rb", "Loading startfile\n", verbose);
	Geometry start_file_geometry;
	NUMERIC_TYPE no_data_value;
	AsciiRaster::read_header(file, start_file_geometry, no_data_value);
	AsciiRaster::match_cell_dimensions_or_die(geometry, start_file_geometry,
			"StartFile::load");
	AsciiRaster::read(file, array, geometry, pitch, offset);
	fclose(file);
	
	zero_no_data_values(array, geometry, pitch, offset, no_data_value);
}

void lis::StartFile::zero_no_data_values
(
	NUMERIC_TYPE* array,
	Geometry& geometry,
	int pitch,
	int offset,
	NUMERIC_TYPE no_data_value
)
{
	for (int j=0; j<geometry.ysz; j++)
	{
		for (int i=0; i<geometry.xsz; i++)
		{
			NUMERIC_TYPE& value = array[j*pitch + i + offset];
			if (FABS(value - no_data_value) < C(1e-5)) value = C(0.0);
		}
	}
}

void lis::StartFile::subtract_dem
(	
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* DEM,
	Geometry& geometry,
	int pitch,
	int offset
)
{
	for (int j=0; j<geometry.ysz; j++)
	{
		for (int i=0; i<geometry.xsz; i++)
		{
			NUMERIC_TYPE& H_value = H[j*pitch + i + offset];
			NUMERIC_TYPE DEM_value = DEM[j*pitch + i + offset];

			H_value = H_value - DEM_value;
		}
	}
}

std::future<void> lis::Snapshot::initialise_async()
{
	return std::async(std::launch::async, []{});
}
