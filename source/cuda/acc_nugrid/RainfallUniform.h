#pragma once

#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../lisflood.h"

namespace lis
{
namespace cuda
{
namespace acc_nugrid
{

typedef struct RainfallUniform
{
	NUMERIC_TYPE        value;
	char                timeseries[32] = { '\0' };
	int                 row = 0;

	RainfallUniform
	(
		const NUMERIC_TYPE& value
	)
	{
		this->value = value;
	}

	void update_rainfall
	(
		const char* rainfilename,
		const NUMERIC_TYPE& time_now
	)
	{
		char str[255];
		char buf[32] = { '\0' };

		FILE* fp = fopen(rainfilename, "r");

		if (NULL == fp)
		{
			fprintf(stderr, "Error opening time varying rainfall %s, file: %s, line: %d.\n", rainfilename, __FILE__, __LINE__);
			exit(-1);
		}

		char* timeseriesptr = timeseries;

		fgets(str, sizeof(str), fp);
		sscanf(str, "%s", buf);

		NUMERIC_TYPE value_1 = C(0.0);
		NUMERIC_TYPE value_2 = C(0.0);

		NUMERIC_TYPE t_1 = C(0.0);
		NUMERIC_TYPE t_2 = C(0.0);

		int num_rows_timeseries = 0;
		int time_multiplier = 1;

		fgets(str, sizeof(str), fp);

		sscanf(str, "%d %s", &num_rows_timeseries, buf);

		time_multiplier = (!strncmp(buf, "seconds", 7)) ? 1 : (!strncmp(buf, "minutes", 7)) ? 60 : 3600;

		for (int i = 0; i < row + 1; i++) fgets(str, sizeof(str), fp);

		sscanf(str, "%" NUM_FMT " %" NUM_FMT, &value_1, &t_1);

		fgets(str, sizeof(str), fp);

		sscanf(str, "%" NUM_FMT " %" NUM_FMT, &value_2, &t_2);

		if (time_now > t_2*time_multiplier)
		{
			row++;

			t_1 = t_2*time_multiplier;
			value_1 = value_2;

			fgets(str, sizeof(str), fp);

			sscanf(str, "%" NUM_FMT " %" NUM_FMT, &value_2, &t_2);
		}

		fclose(fp);

		if (row + 1 >= num_rows_timeseries)
		{
			row--;
			value = C(0.0);
		}
		else {
			value = value_1 + (value_2 - value_1) / ((t_2 - t_1) * time_multiplier) * (time_now - t_1 * time_multiplier);
			value /= (C(1000.0) * C(3600.0));
		}
	}

} RainfallUniform;

}
}
}