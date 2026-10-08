#pragma once
#include "../lisflood.h"
#include <cstdio>

namespace lis
{

typedef struct MassStats
{
	NUMERIC_TYPE in;
	NUMERIC_TYPE out;
} MassStats;

typedef struct StatsEntry
{
	NUMERIC_TYPE t;
	NUMERIC_TYPE dt;
	long iteration;
	NUMERIC_TYPE min_dt;
	NUMERIC_TYPE area;
	NUMERIC_TYPE volume;
	MassStats mass;
	NUMERIC_TYPE volume_error;
	NUMERIC_TYPE discharge_error;
} StatsEntry;

class Stats
{
public:
	Stats
	(
		const char* resroot,
		NUMERIC_TYPE interval,
		NUMERIC_TYPE next_save,
		int checkpoint,
		NUMERIC_TYPE time
	);

	void write_header(NUMERIC_TYPE t);

	bool need_to_write
	(
		NUMERIC_TYPE t
	);

	void write
	(
		StatsEntry entry
	);

	NUMERIC_TYPE next_save;

	~Stats();

private:
	FILE* file;
	NUMERIC_TYPE interval;
	
};

}
