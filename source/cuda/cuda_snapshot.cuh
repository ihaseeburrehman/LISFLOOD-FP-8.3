#pragma once
#include "../lisflood.h"
#include "../geometry.h"
#include "io.h"
#include <future>

namespace lis
{
namespace cuda
{

template<typename F>
class Snapshot
{
public:
	Snapshot
	(
		const char* resrootname,
		NUMERIC_TYPE interval,
		NUMERIC_TYPE next_save,
		int counter,
		F& U,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose,
		int precision = DEFAULT_PRECISION
	);

	void write_if_needed
	(
		F& d_U,
		NUMERIC_TYPE t,
		const int& call_gzip
	);

	void write_max_files
	(
		F& d_U
	);

	void enable_elevation_writer();

	void enable_discharge_writer();

	void enable_velocity_writer
	(
		NUMERIC_TYPE DepthThresh
	);

	virtual void wait() = 0;

	int counter;
	NUMERIC_TYPE next_save;

protected:
	virtual void write(const int& call_gzip) = 0;

	virtual void write_maxes() = 0;

	template<typename T>
	std::future<void> write_async
	(
		T field,
		const char* suffix,
		int outflag,
		const int& call_gzip
	);

	template<typename T>
	void write_max
	(
		T field,
		const char* suffix,
		int outflag
	);

	F& U;
	bool write_elevation = false;
	bool write_discharge = false;
	bool write_velocity = false;
	NUMERIC_TYPE DepthThresh = C(0.0);

private:
	const char* resrootname;
	Geometry& geometry;
	NUMERIC_TYPE interval;
	
	int pitch;
	int offset;
	int verbose;
	int precision;

};

class VelocityWriter
{
public:
	VelocityWriter
	(
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* discharge,
		NUMERIC_TYPE DepthThresh
	);

	NUMERIC_TYPE operator[]
	(
		int idx
	);

private:
	NUMERIC_TYPE* H;
	NUMERIC_TYPE* discharge;
	NUMERIC_TYPE DepthThresh;
};

class ElevationWriter
{
public:
	ElevationWriter
	(
		const NUMERIC_TYPE* H,
		const NUMERIC_TYPE* DEM
	);

	NUMERIC_TYPE operator[]
	(
		int idx
	) const;

private:
	const NUMERIC_TYPE* H;
	const NUMERIC_TYPE* DEM;
};

}
}

#include "cuda_snapshot.templates.cu"
