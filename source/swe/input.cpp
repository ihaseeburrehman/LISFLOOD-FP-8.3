#include "input.h"

void read_ascfile
(
	const char* filename,
	Pars *Parptr,
	NUMERIC_TYPE *field,
	const char* message,
	const int verbose
)
{
	FILE *fp;
	char dummy[800];
	NUMERIC_TYPE no_data_value = -9999;

	fp = fopen_or_die(filename, "r", message, verbose);
	for (int i=0; i<5; i++)
	{
		fscanf(fp, "%s %s", dummy, dummy);
	}
	fscanf(fp, "%s %" NUM_FMT"", dummy, &no_data_value);
	for (int j=0; j<Parptr->ysz; j++)
	{
		for (int i=0; i<Parptr->xsz; i++)
		{
			fscanf(fp, "%" NUM_FMT"", field + j*Parptr->xsz + i);
			if (FABS(field[i + j*Parptr->xsz] - no_data_value) < C(1e-12))
			{
				field[j*Parptr->xsz + i] = C(0.0);
			}
		}
	}
	fclose(fp);
}
