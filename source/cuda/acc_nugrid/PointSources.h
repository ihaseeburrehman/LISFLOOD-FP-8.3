#pragma once

#include "cuda_utils.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MortonCode.h"
#include "InletTypes.h"
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct SourceAndRow
{
	NUMERIC_TYPE src;
	int  row;

} SourceAndRow;

typedef struct PointSources
{
	MortonCode* h_codes;
	MortonCode* d_codes;
	NUMERIC_TYPE*       h_srcs;
	NUMERIC_TYPE*       d_srcs;
	int         num_srcs;
	int*        h_src_types;
	int*        d_src_types;
	char*       timeseries;
	int*        rows;
	bool        is_copy = false;

	PointSources(const int& num_srcs)
	{
		size_t bytes_codes = sizeof(MortonCode) * num_srcs;
		size_t bytes_srcs = sizeof(NUMERIC_TYPE) * num_srcs;
		size_t bytes_src_types = sizeof(int) * num_srcs;

		h_codes = (num_srcs > 0) ? new MortonCode[num_srcs] : nullptr;
		d_codes = (num_srcs > 0) ? (MortonCode*)malloc_device(bytes_codes) : nullptr;
		h_srcs = (num_srcs > 0) ? new NUMERIC_TYPE[num_srcs]() : nullptr;
		d_srcs = (num_srcs > 0) ? (NUMERIC_TYPE*)malloc_device(bytes_srcs) : nullptr;
		h_src_types = (num_srcs > 0) ? new int[num_srcs] : nullptr;
		d_src_types = (num_srcs > 0) ? (int*)malloc_device(bytes_src_types) : nullptr;
		timeseries = (num_srcs > 0) ? new char[num_srcs * 32]() : nullptr;
		rows = (num_srcs > 0) ? new int[num_srcs]() : nullptr;

		this->num_srcs = num_srcs;
	}

	PointSources(const PointSources& original) { *this = original; is_copy = true; }

	~PointSources()
	{
		if (!is_copy)
		{
			if (h_codes != nullptr) delete[] h_codes;
			if (h_srcs != nullptr) delete[] h_srcs;
			if (h_src_types != nullptr) delete[] h_src_types;

			if (h_codes != nullptr) free_device(d_codes);
			if (d_srcs != nullptr) free_device(d_srcs);
			if (d_src_types != nullptr) free_device(d_src_types);
		}
	}

	SourceAndRow update_source //TOREMOVE
	(
		const char* srcfilename,
		const NUMERIC_TYPE& time_now,
		const int&  row_in,
		const int&  src_type,
		const char* timeseries
	)
	{
		if (src_type != HVAR && src_type != QVAR) return { C(0.0), 0 };
		
		char str[255];
		char buf[32] = {'\0'};
		
		FILE* fp = fopen(srcfilename, "r");

		if (NULL == fp)
		{
			fprintf(stderr, "Error opening time varying boundary condition file %s, file: %s, line: %d.\n", srcfilename, __FILE__, __LINE__);
			exit(-1);
		}

		int num_char_timeseries = 0;

		while ( *(timeseries + num_char_timeseries) != '\0' ) num_char_timeseries++;

		while ( strncmp(buf, timeseries, num_char_timeseries) )
		{
			if ( NULL == fgets(str, sizeof(str), fp) )
			{
				fprintf(stderr, "Error reading point source time series \"%s\", file: %s, line: %d.\n", timeseries, __FILE__, __LINE__);
				fclose(fp);
				exit(-1);
			}

			sscanf(str, "%s", buf);
		}
	
		NUMERIC_TYPE src   = C(0.0);
		NUMERIC_TYPE src_1 = C(0.0);
		NUMERIC_TYPE src_2 = C(0.0);

		NUMERIC_TYPE t_1 = C(0.0);
		NUMERIC_TYPE t_2 = C(0.0);

		int row                 = row_in;
		int num_rows_timeseries = 0;
		int time_multiplier     = 1;

		fgets(str, sizeof(str), fp);

		sscanf(str, "%d %s", &num_rows_timeseries, buf);

		time_multiplier = ( !strncmp(buf, "seconds", 7) ) ? 1 : ( strncmp(buf, "minutes", 7) ) ? 60 : 3600;

		for (int i = 0; i < row + 1; i++) fgets(str, sizeof(str), fp);

		sscanf(str, "%" NUM_FMT " %" NUM_FMT, &src_1, &t_1);

		fgets(str, sizeof(str), fp);

		sscanf(str, "%" NUM_FMT " %" NUM_FMT, &src_2, &t_2);

		if (time_now > t_2)
		{
			row++;

			if (row + 1 == num_rows_timeseries)
			{
				fprintf
				(
					stderr,
					"Error: simulation time exceeded time specified in time series \"%s\", file: %s, line: %d.\n", 
					timeseries, __FILE__, __LINE__
				);

				fclose(fp);
				exit(-1);
			}

			t_1   = t_2;
			src_1 = src_2;

			fgets(str, sizeof(str), fp);

			sscanf(str, "%" NUM_FMT " %" NUM_FMT, &src_2, &t_2);
		}

		fclose(fp);

		src = src_1 + (src_2 - src_1) / ( (t_2 - t_1) * time_multiplier ) * (time_now - t_1);

		return { src, row };
	}

	void update_all_sources
	(
		const char* input_filename, 
		const NUMERIC_TYPE& time_now,
		const BoundCs& h_src
	)
	{
		for (int i = 0; i < this->num_srcs; i++)
		{
			//int   row_in     =  this->rows[i];
			int   src_type   =  this->h_src_types[i];
			//char* timeseries = &this->timeseries[i * 32];

			//SourceAndRow src_and_row = this->update_source
			//(
			//	input_filename, 
			//	time_now, 
			//	row_in, 
			//	src_type, 
			//	timeseries
			//);

			if (src_type == HVAR || src_type == QVAR)
			{
//				this->h_srcs[i] = src_and_row.src;
//				this->rows[i]   = src_and_row.row;
				this->h_srcs[i] = InterpolateTimeSeries(h_src.PS_TimeSeries[i], time_now);

			}
		}

//		size_t bytes_srcs = sizeof(NUMERIC_TYPE) * this->num_srcs;

//		copy_cuda(this->d_srcs, this->h_srcs, bytes_srcs);
		cudaError_t error = cudaMemcpy(d_srcs, h_srcs, sizeof(NUMERIC_TYPE) * this->num_srcs, cudaMemcpyDefault); // removed teh wrapper
	}

	__device__ __forceinline__
	NUMERIC_TYPE q_src(const NUMERIC_TYPE& dt, const NUMERIC_TYPE& dx, const int idx) { return this->d_srcs[idx] * dt / dx; }

} PointSources;

}
}
}