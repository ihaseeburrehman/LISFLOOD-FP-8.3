#pragma once
#define DEFAULT_SAMPLE_BUFFER_SIZE 1
#include "../lisflood.h"
#include "sample.h"

namespace lis
{
namespace cuda
{

class Sampler
{
public:
	Sampler
	(
		SamplePoints& d_sample_points,
		int& sample_buf_idx,
		int verbose
	);

	void load_sample_points
	(
		const char* filename,
		Geometry& geometry,
		int pitch,
		int offset
	);

	void open_stage_file
	(
	 	const char* filename,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void open_gauge_file
	(
	 	const char* filename,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void write_stage_header
	(
		NUMERIC_TYPE* DEM,
		const char* sample_points_filename,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void write_gauge_header
	(
		NUMERIC_TYPE* DEM,
		const char* sample_points_filename,
		NUMERIC_TYPE t
	);

	void sample
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* HU,
		NUMERIC_TYPE* HV,
		NUMERIC_TYPE t
	);
	
	void sample_ACC
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* Vx,
		NUMERIC_TYPE* Vy,
		NUMERIC_TYPE t
	);

	void write_if_buffer_full();

	void write();
	
	~Sampler();

private:
	bool buffer_full();

	void copy_buffer();

	void initialise_sample_points();

	void free
	(
		SamplePoints& d_sample_points
	);

	void allocate_pinned
	(
	 	SampleBuffer& buf,
		int points,
		int size = DEFAULT_SAMPLE_BUFFER_SIZE
	);

	void allocate_device
	(
	 	SampleBuffer& buf,
		int points,
		int size = DEFAULT_SAMPLE_BUFFER_SIZE
	);

	void free_pinned
	(
		SampleBuffer& buf
	);

	void free_device
	(
		SampleBuffer& buf
	);

	bool active = false;
	bool write_speed = false;
	int& sample_buf_idx;
	SamplePoints sample_points;
	SamplePoints& d_sample_points;
	SampleBuffer sample_buf;
	SampleBuffer d_sample_buf;
	StageFile stage_file;
	GaugeFile gauge_file;
	int verbose;
};

}
}
