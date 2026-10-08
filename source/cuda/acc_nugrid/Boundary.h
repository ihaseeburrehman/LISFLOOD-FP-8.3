#pragma once

#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Directions.h"
#include "MortonCode.h"
#include "InletTypes.h"
#include "../../lisflood.h"
#include "generate_morton_code.cuh"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct Boundary
{
	int  start          = 0;
	int  end            = 0;
	MortonCode* codes          = nullptr;
	int         bdytype        = CLOSED;
	NUMERIC_TYPE        inlet          = C(0.0);
	char        timeseries[32] = {'\0'};
	int         timeseries_len = 0;
	NUMERIC_TYPE*       time_data      = nullptr;
	NUMERIC_TYPE*       inlet_data     = nullptr;
	int         row            = 0;
	int         direction;
	bool        is_copy        = false;

	Boundary
	(
		const Fnames& filenames,
		const Pars& pars,
		const int                   direction
	)
	:
		direction(direction)
	{
		const int num_cells = read_num_cells
		(
			filenames,
			pars
		);
		
		codes = (num_cells > 0) ? new MortonCode[num_cells] : nullptr;

		read_bdy_conds
		(
			filenames,
			pars
		);
		
		read_time_series(filenames);
	}

	Boundary(const Boundary& original) { *this = original; is_copy = true; }

	~Boundary()
	{
		if (!is_copy)
		{
			if (codes      != nullptr) delete[] codes;
			if (time_data  != nullptr) delete[] time_data;
			if (inlet_data != nullptr) delete[] inlet_data;
		}
	}

	__device__ __forceinline__
	bool bound(const int& coordinate) const { return (coordinate >= start && coordinate <= end); }

	__host__ __forceinline__
	int num_cells() const { return end - start + 1; }

	__device__ __forceinline__
		NUMERIC_TYPE q_src(const NUMERIC_TYPE& dt, const NUMERIC_TYPE& dx) { return inlet * dt / dx; }

	int read_num_cells
	(
		const Fnames& filenames,
		const Pars& pars
	)
	{
		
		NUMERIC_TYPE upper = C(0.0);
		NUMERIC_TYPE lower = C(0.0);

		char bcidir  = '\0';
		char filedir = '\0';

		char str[255];
		char bdytype_buf[8] = {'\0'};

		NUMERIC_TYPE origin = C(0.0);

		switch (direction)
		{
		case NORTH:
			origin = pars.blx;
			bcidir = 'N';
			break;
		case EAST:
			origin = pars.bly;
			bcidir = 'E';
			break;
		case SOUTH:
			origin = pars.blx;
			bcidir = 'S';
			break;
		case WEST:
			origin = pars.bly;
			bcidir = 'W';
			break;
		default:
			break;
		}
	

		if (strlen(filenames.bcifilename) == 0)
		{
			printf("No enforced boundary cells counted for boundary.\n");
			return 0;
		}

		FILE* fp = fopen(filenames.bcifilename, "r");

		if (NULL == fp)
		{
			fprintf(stderr, "Error opening boundary condition file, file: %s, line: %d.\n", __FILE__, __LINE__);
			exit(-1);
		}

		while (filedir != bcidir)
		{

			if ( NULL == fgets(str, sizeof(str), fp) )
			{
				fprintf(stdout, "No enforced boundary cells counted for boundary %c.\n", bcidir);
				fclose(fp);

				return 0;
			}

			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s", &filedir, &lower, &upper, bdytype_buf);
		}

		fclose(fp);

		int start = (lower - origin) / pars.dx;
		int end = (upper - origin) / pars.dx - 1;

		return end - start + 1;
	}

	void gen_bdy_morton_codes
	(
		const Pars& pars
	)
	{
		int current = 0;
	
		switch (direction)
		{
		case SOUTH:
		{
			for (int i = 0; i < this->num_cells(); i++)
			{
				current = this->start + i;

				this->codes[i] = generate_morton_code(current, 0);
			}

			break;
		}
		case NORTH:
		{
			for (int i = 0; i < this->num_cells(); i++)
			{
				current = this->start + i;

				this->codes[i] = generate_morton_code(current, pars.ysz - 1);
			}

			break;
		}
		case EAST:
		{
			for (int i = 0; i < this->num_cells(); i++)
			{
				current = this->start + i;

				this->codes[i] = generate_morton_code(pars.xsz - 1, current);
			}

			break;
		}
		case WEST:
		{
			for (int i = 0; i < this->num_cells(); i++)
			{
				current = this->start + i;

				this->codes[i] = generate_morton_code(0, current);
			}

			break;
		}
		default:
			break;
		}
	}

	void read_bdy_conds
	(
		const Fnames& filenames,
		const Pars& pars
	)
	{

		NUMERIC_TYPE upper = C(0.0);
		NUMERIC_TYPE lower = C(0.0);

		char bcidir  = '\0';
		char filedir = '\0';

		char str[255];
		char bdytype_buf[8] = {'\0'};

		NUMERIC_TYPE origin = C(0.0);

		switch (direction)
		{
		case NORTH:
			origin = pars.blx;
			bcidir = 'N';
			break;
		case EAST:
			origin = pars.bly;
			bcidir = 'E';
			break;
		case SOUTH:
			origin = pars.blx;
			bcidir = 'S';
			break;
		case WEST:
			origin = pars.bly;
			bcidir = 'W';
			break;
		default:
			break;
		}
	
		if (strlen(filenames.bcifilename) == 0)
		{
			printf("No specifications found for boundary, proceeding with closed boundaries.\n");
			return;
		}

		FILE* fp = fopen(filenames.bcifilename, "r");

		if (NULL == fp)
		{
			fprintf(stderr, "Error opening boundary condition file, file: %s, line: %d.\n", __FILE__, __LINE__);
			exit(-1);
		}
		
		int  bdytype        = CLOSED;
		NUMERIC_TYPE inlet          = C(0.0);
		char timeseries[32] = {'\0'};

		while (filedir != bcidir)
		{
			if ( NULL == fgets(str, sizeof(str), fp) )
			{
				fprintf(stdout, "No specifications found for boundary %c, proceeding with closed boundaries.\n", bcidir);
				fclose(fp);

				return;
			}

			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s", &filedir, &lower, &upper, bdytype_buf);
		}

		fclose(fp);
	
		if (bcidir == 'N') {
//			int BCi1 = -1;

			if (lower < pars.blx) lower = pars.blx;
			if (lower > pars.blx + pars.xsz * pars.dx) lower = pars.blx + pars.xsz * pars.dx;
			if (upper < pars.blx) upper = pars.blx;
			if (upper > pars.blx + pars.xsz * pars.dx) upper = pars.blx + pars.xsz * pars.dx;

		}

		if (bcidir == 'W') {

//			int BCi1 = -1;

			if (lower < pars.bly) lower = pars.bly;
			if (lower > pars.tly) lower = pars.tly;
			if (upper < pars.bly) upper = pars.bly;
			if (upper > pars.tly) upper = pars.tly;

		}

		if (bcidir == 'S')
		{
//			int BCi1 = -1;

			if (lower < pars.blx) lower = pars.blx;
			if (lower > pars.blx + pars.xsz * pars.dx) lower = pars.blx + pars.xsz * pars.dx;
			if (upper < pars.blx) upper = pars.blx;
			if (upper > pars.blx + pars.xsz * pars.dx) upper = pars.blx + pars.xsz * pars.dx;

		}

		if (bcidir == 'E')
		{
//			int BCi1 = -1;

			if (lower < pars.bly) lower = pars.bly;
			if (lower > pars.tly) lower = pars.tly;
			if (upper < pars.bly) upper = pars.bly;
			if (upper > pars.tly) upper = pars.tly;

		}


		if ( !strncmp(bdytype_buf, "CLOSED", 6) )
		{
			bdytype = CLOSED;
		}
		else if ( !strncmp(bdytype_buf, "FREE", 4) )
		{
			bdytype = FREE;
			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &filedir, &lower, &upper, bdytype_buf, &inlet);
		}
		else if ( !strncmp(bdytype_buf, "HFIX", 4) )
		{
			bdytype = HFIX;
			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &filedir, &lower, &upper, bdytype_buf, &inlet);
		}
		else if ( !strncmp(bdytype_buf, "HVAR", 4) )
		{
			bdytype = HVAR;
			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &filedir, &lower, &upper, bdytype_buf, timeseries);
		}
		else if ( !strncmp(bdytype_buf, "QFIX", 4) )
		{
			bdytype = QFIX;
			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &filedir, &lower, &upper, bdytype_buf, &inlet);
		}
		else if ( !strncmp(bdytype_buf, "QVAR", 4) )
		{
			bdytype = QVAR;
			sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &filedir, &lower, &upper, bdytype_buf, timeseries);
		}

		this->start   = (lower - origin) / pars.dx;
		this->end     = (upper - origin) / pars.dx - 1;
		gen_bdy_morton_codes(pars);
		this->bdytype = bdytype;
		this->inlet   = inlet;
		sprintf(this->timeseries, "%s", timeseries);
	}

	void read_time_series
	(
		const Fnames& filenames
	)
	{
		if ( (bdytype != HVAR && bdytype != QVAR) )
		{
			time_data  = nullptr;
			inlet_data = nullptr;
			
			return;
		}

		char str[255]        = {'\0'};
		char buf[64]         = {'\0'};

		FILE* fp = fopen(filenames.bdyfilename, "r");

		if (NULL == fp)
		{
			fprintf(stderr, "Error opening time varying boundary condition file %s, file: %s, line: %d.\n", filenames.bdyfilename, __FILE__, __LINE__);
			exit(-1);
		}

		char* timeseriesptr = timeseries;

		int num_char_timeseries = 0;

		while (*(timeseriesptr + num_char_timeseries) != '\0') num_char_timeseries++;

		while ( strncmp(buf, timeseries, num_char_timeseries) )
		{
			if ( NULL == fgets(str, sizeof(str), fp) )
			{
				fprintf(stderr, "Error reading boundary inlet time series \"%s\", file: %s, line: %d.\n", timeseries, __FILE__, __LINE__);
				fclose(fp);
				exit(-1);
			}

			sscanf(str, "%s", buf);
		}

		int num_rows_timeseries = 0;

		fgets(str, sizeof(str), fp);

		sscanf(str, "%d %s", &num_rows_timeseries, buf);

		int time_multiplier = ( !strncmp(buf, "seconds", 7) ) ? 1 : ( !strncmp(buf, "minutes", 7) ) ? 60 : 3600;

		timeseries_len = num_rows_timeseries;

		if (num_rows_timeseries == 0)
		{
			fprintf(stderr, "Zero entries for timeseries: \"%s\", file: %s, line: %d.\n", timeseries, __FILE__, __LINE__);
			fclose(fp);
			exit(-1);
		}

		time_data  = (num_rows_timeseries > 0) ? new NUMERIC_TYPE[num_rows_timeseries] : nullptr;
		inlet_data = (num_rows_timeseries > 0) ? new NUMERIC_TYPE[num_rows_timeseries] : nullptr;

		for (int i = 0; i < num_rows_timeseries; i++)
		{
			fgets(str, sizeof(str), fp);

			sscanf(str, "%" NUM_FMT " %" NUM_FMT, &inlet_data[i], &time_data[i]);

			time_data[i] *= time_multiplier;
		}

		fclose(fp);
	}

	void update_inlet
	(
		const NUMERIC_TYPE& time_now
	)
	{
		if ( (bdytype != HVAR && bdytype != QVAR) ) return;

		if ( (row - 1) < timeseries_len )
		{
			NUMERIC_TYPE t_1 = time_data[row];
			NUMERIC_TYPE t_2 = time_data[row + 1];

			if (time_now > t_2)
			{
				row++; 
				
				t_1 = time_data[row];
				t_2 = time_data[row + 1];
			}

			NUMERIC_TYPE inlet_1 = inlet_data[row];
			NUMERIC_TYPE inlet_2 = inlet_data[row + 1];

			inlet = inlet_1 + (inlet_2 - inlet_1) / (t_2 - t_1) * (time_now - t_1);
		}
		else
		{
			inlet = inlet_data[timeseries_len - 1];
		}
	}

} Boundary;

}
}
}