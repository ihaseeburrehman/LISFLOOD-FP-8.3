#pragma once
#include "../lisflood.h"
#include "../geometry.h"
#include "../rain/rain.h"
//#include "acc/cuda_acc_flow.cuh"
//#include "cuda_max_field.cuh"
#include "stats.h"
#include "cuda_stats.cuh"

namespace lis
{
namespace cuda
{

class Simulation
{
protected:
	void print_device_info();

	void initialise_H
	(
		NUMERIC_TYPE* H,
		const char* filename,
		States& states,
		NUMERIC_TYPE* DEM,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void nullify_max
	(
		NUMERIC_TYPE* maxHtm,
		NUMERIC_TYPE* initHtm,
		NUMERIC_TYPE* totalHtm,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void initialise_discharge
	(
		NUMERIC_TYPE* HU,
		NUMERIC_TYPE* HV,
		const char* filename,
		States& states,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void initialise_manning
	(
		NUMERIC_TYPE* manning,
		const char* filename,
		Geometry& geometry,
		int pitch,
		int offset,
		int verbose
	);

	void update_geometry
	(
		Pars& dst,
		Geometry& src
	);

	void load_boundaries
	(
		Fnames& filenames,
		States& states,
		Pars& pars,
		BoundCs& boundCs,
		int verbose
	);

	void read_checkpoint
	(
		Fnames& filenames,
		States& states,
		Pars& pars,
		Solver& solver,
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* Qx,
		NUMERIC_TYPE* Qy,
		NUMERIC_TYPE*  maxH,
		NUMERIC_TYPE* totalHtm, 
		NUMERIC_TYPE*  maxHtm, 
		NUMERIC_TYPE* initHtm,
		StatsCollector& Stats_Entry,
		int verbose
	);

	void write_checkpoint
	(
		Fnames& filenames,
		States& states,
		Pars& pars,
		Solver& solver,
		NUMERIC_TYPE* H,
		NUMERIC_TYPE* Qx,
		NUMERIC_TYPE* Qy,
		NUMERIC_TYPE* maxH,
		NUMERIC_TYPE* totalHtm,
		NUMERIC_TYPE* maxHtm,
		NUMERIC_TYPE* initHtm,
		lis::StatsEntry& Stats_Entry,
		NUMERIC_TYPE dt,
		int verbose
	);


};

}
}
