#pragma once
#include "../geometry.h"

namespace lis
{

typedef struct SamplePoints
{
	int count;
	NUMERIC_TYPE* x;
	NUMERIC_TYPE* y;
	int* idx;
	int* idx_Qx1; 
	int* idx_Qx2; 
	int* idx_Qy1; 
	int* idx_Qy2; 
	bool* inside_domain;
} SamplePoints;

typedef struct SampleBuffer
{
	int size;
	NUMERIC_TYPE* time;
	NUMERIC_TYPE* H;
	NUMERIC_TYPE* speed;
} SampleBuffer;

struct Sample
{
	static void initialise
	(
		SamplePoints& sample_points,
		const char* filename,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	static void free
	(
		SamplePoints& sample_points
	);
};

class StageFile
{
public:
	StageFile
	(
		SamplePoints& sample_points
	);

	void open
	(
		const char* resroot,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void write_header
	(
		NUMERIC_TYPE* DEM,
		const char* sample_points_filename,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void write
	(
		SampleBuffer& sample_buf,
		int sample_buf_idx
	);

	~StageFile();

private:
	FILE* file = nullptr;
	const SamplePoints& sample_points;
};

class GaugeFile
{
public:
	GaugeFile
	(
		SamplePoints& sample_points
	);

	void open
	(
		const char* resroot,
		int checkpoint,
		NUMERIC_TYPE t
	);

	void write_header
	(
		NUMERIC_TYPE* DEM,
		const char* sample_points_filename,
		NUMERIC_TYPE t
	);

	void write
	(
		SampleBuffer& sample_buf,
		int sample_buf_idx
	);

	~GaugeFile();

private:
	FILE* file = nullptr;
	const SamplePoints& sample_points;
};

}
