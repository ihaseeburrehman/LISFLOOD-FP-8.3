#pragma once
#include <functional>
#include "../../lisflood.h"

namespace dg2
{
	void initialise_field
	(
		std::function<NUMERIC_TYPE(NUMERIC_TYPE, NUMERIC_TYPE)> func,
		Pars *Parptr,
		NUMERIC_TYPE *F,
		NUMERIC_TYPE *F1x,
		NUMERIC_TYPE *F1y
	);

	void initialise_field
	(
		std::function<NUMERIC_TYPE(NUMERIC_TYPE, NUMERIC_TYPE)> func,
		NUMERIC_TYPE no_data,
		Pars *Parptr,
		NUMERIC_TYPE *F,
		NUMERIC_TYPE *F1x,
		NUMERIC_TYPE *F1y
	);

	void downscale
	(
		NUMERIC_TYPE* wd0,
		NUMERIC_TYPE* wd1x,
		NUMERIC_TYPE* wd1y,
		NUMERIC_TYPE no_data,
		Pars* Parptr,
		NUMERIC_TYPE* F
	);

	void initialise_h_from_eta
	(
		NUMERIC_TYPE eta,
		Pars *Parptr,
		Arrays *Arrptr
	);

	void read_dem_slopes
	(
		Fnames *Fnameptr,
		Pars *Parptr,
		Arrays *Arrptr,
		const int verbose
	);

    void zero_dem_perimeter_slopes
    (
		Pars *Parptr,
		Arrays *Arrptr
    );

	void read_h_slopes
	(
		Fnames *Fnameptr,
		Pars *Parptr,
		Arrays *Arrptr,
		const int verbose
	);

	void read_discharge_slopes
	(
		Fnames *Fnameptr,
		Pars *Parptr,
		Arrays *Arrptr,
		const int verbose
	);

	void write_dem
	(
		const char *root,
		Pars *Parptr,
		Arrays *Arrptr
	);

	void write_DScaled
	(
		const char* root,
		Pars* Parptr,
		Arrays* Arrptr
	);

	void write_startfile
	(
		const char *root,
		Pars *Parptr,
		Arrays *Arrptr
	);

	void allocate_fields
	(
		Pars *Parptr,
		Arrays *Arrptr
	);

	void allocate
	(
		BoundaryValues& boundary,
		Pars *Parptr
	);

	void deallocate_fields
    (
		Pars *Parptr,
        Arrays *Arrptr
    );

	void deallocate(BoundaryValues& boundary);
}
