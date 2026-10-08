#include <cstdlib>
#include "../swe/dg2/fields.h"
#include "../lisflood.h"
#include "../utility.h"

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
	ReadConfiguration(argc, argv, Fnameptr, Statesptr, Parptr, Solverptr,
			verbose);

	char startfilename[256];
	strcpy(startfilename, Fnameptr->startfilename);
	strcat(startfilename, ".raw");

	NUMERIC_TYPE no_data_value = NULLVAL;
	FILE *fp = LoadDomainGeometry(startfilename, Parptr, verbose, no_data_value);
	NUMERIC_TYPE *rawdem;
	rawdem = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	LoadDEMData(Parptr, rawdem, fp, no_data_value);
	fclose(fp);

	Arrays arrays;
	Arrays *Arrptr = &arrays;

	Arrptr->DEM = memory_allocate_numeric_legacy(Parptr->xsz*Parptr->ysz);
	dg2::allocate_fields(Parptr, Arrptr);

	DG2DEM converter(Parptr, rawdem);

	dg2::initialise_field(converter, Parptr->nodata_elevation, Parptr,
			Arrptr->DEM, Arrptr->DEM1x, Arrptr->DEM1y);

	dg2::write_dem(Fnameptr->startfilename, Parptr, Arrptr);
	
	memory_free_legacy(&Arrptr->DEM);
	dg2::deallocate_fields(Parptr, Arrptr);
	return EXIT_SUCCESS;
}

