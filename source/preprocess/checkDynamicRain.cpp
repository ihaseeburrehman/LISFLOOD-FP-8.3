#include <cstdlib>
#include "../lisflood.h"
#include "../utility.h"
#include "../rain/rain.h"

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
    
	NUMERIC_TYPE no_data_value = NULLVAL;
	FILE *fp = LoadDomainGeometry(Fnameptr->demfilename, Parptr, verbose,
            no_data_value);
	NUMERIC_TYPE *DEM;
	DEM = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	LoadDEMData(Parptr, DEM, fp, no_data_value);
	fclose(fp);

    DynamicRain rain(fnames.dynamicrainfilename, verbose);
    printf("has_same_origin: ");
    printf(rain.has_same_origin(Parptr) ? "true" : "false");
    printf("\n");
    printf("is_tile_size_multiple_of_grid: ");
    printf(rain.is_tile_size_multiple_of_grid(Parptr) ? "true" : "false");
    printf("\n");

    return EXIT_SUCCESS;
}
