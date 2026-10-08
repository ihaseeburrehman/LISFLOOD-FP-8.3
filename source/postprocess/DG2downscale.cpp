#include <cstdlib>
#include "../swe/dg2/fields.h"
#include "../lisflood.h"
#include "../utility.h"
#include <string>

class DG2DEM
{
private:
	Pars *Parptr;
	NUMERIC_TYPE *rawdem;

public:
	DG2DEM
	(
		Pars *Parptr,
		NUMERIC_TYPE *rawdem
	)
	:
	Parptr(Parptr),
	rawdem(rawdem)
	{}

	NUMERIC_TYPE operator()(NUMERIC_TYPE x, NUMERIC_TYPE y)
	{
		int i = static_cast<int>((x - Parptr->blx)/Parptr->dx);
		int j = static_cast<int>((Parptr->tly - y)/Parptr->dy);

		i = FMIN(i, Parptr->xsz-1);
		j = FMIN(j, Parptr->ysz-1);

		return rawdem[j*Parptr->xsz + i];
	}
};

int main(int argc, char *argv[])
{
	Fnames fnames;
	Fnames *Fnameptr = &fnames;

	States states;
	States *Statesptr = &states;

	Pars pars;
	Pars *Parptr = &pars;

	Solver solver;
	Solver *Solverptr = &solver;

	int verbose = ReadVerboseMode(argc, argv);

	char wd0filename[256];
	char wd1xfilename[256];
	char wd1yfilename[256];

	strcpy(wd0filename, argv[argc - 1]);
	strcpy(wd1xfilename, argv[argc - 1]);
	strcpy(wd1yfilename, argv[argc - 1]);
	strcat(wd1xfilename, "1x");
	strcat(wd1yfilename, "1y");


		NUMERIC_TYPE no_data_value = NULLVAL;

		FILE* fp = LoadDomainGeometry(wd0filename, Parptr, verbose, no_data_value);
		FILE* fp1x = LoadDomainGeometry(wd1xfilename, Parptr, verbose, no_data_value);
		FILE* fp1y = LoadDomainGeometry(wd1yfilename, Parptr, verbose, no_data_value);

		NUMERIC_TYPE* wd0;
		NUMERIC_TYPE* wd1x;
		NUMERIC_TYPE* wd1y;


		wd0 = memory_allocate_zero_numeric_legacy(Parptr->xsz * Parptr->ysz);
		wd1x = memory_allocate_zero_numeric_legacy(Parptr->xsz * Parptr->ysz);
		wd1y = memory_allocate_zero_numeric_legacy(Parptr->xsz * Parptr->ysz);


		LoadDEMData(Parptr, wd0, fp, no_data_value);
		LoadDEMData(Parptr, wd1x, fp1x, no_data_value);
		LoadDEMData(Parptr, wd1y, fp1y, no_data_value);
		fclose(fp);
		fclose(fp1x);
		fclose(fp1y);

		Arrays arrays;
		Arrays* Arrptr = &arrays;

		Arrptr->H = memory_allocate_numeric_legacy(4 * Parptr->xsz * Parptr->ysz);


		dg2::downscale(wd0, wd1x, wd1y, Parptr->nodata_elevation, Parptr,
			Arrptr->H);


		dg2::write_DScaled(wd0filename, Parptr, Arrptr);

		memory_free_legacy(&Arrptr->H);

	return EXIT_SUCCESS;
}

