#include "fields.h"
#include "../input.h"
#include "../fields.h"
#include "../../utility.h"

void dg2::initialise_field
(
	std::function<NUMERIC_TYPE(NUMERIC_TYPE, NUMERIC_TYPE)> func,
	Pars *Parptr,
	NUMERIC_TYPE *F,
	NUMERIC_TYPE *F1x,
	NUMERIC_TYPE *F1y
)
{
	initialise_field(func, C(-9999.0), Parptr, F, F1x, F1y);
}

void dg2::initialise_field
(
	std::function<NUMERIC_TYPE(NUMERIC_TYPE, NUMERIC_TYPE)> func,
	NUMERIC_TYPE no_data,
	Pars *Parptr,
	NUMERIC_TYPE *F,
	NUMERIC_TYPE *F1x,
	NUMERIC_TYPE *F1y
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE NW = func(x_vertex(Parptr, i), y_vertex(Parptr, j));
			NUMERIC_TYPE NE = func(x_vertex(Parptr, i+1), y_vertex(Parptr, j));
			NUMERIC_TYPE SW = func(x_vertex(Parptr, i), y_vertex(Parptr, j+1));
			NUMERIC_TYPE SE = func(x_vertex(Parptr, i+1), y_vertex(Parptr, j+1));

			if (FABS(NW - no_data) < 1e-6 || FABS(NE - no_data) < 1e-6
					|| FABS(SW - no_data) < 1e-6 || FABS(SE - no_data) < 1e-6)
			{
				F[j*Parptr->xsz + i] = NULLVAL;
				F1x[j*Parptr->xsz + i] = C(0.0);
				F1y[j*Parptr->xsz + i] = C(0.0);
			}
			else
			{
				F[j*Parptr->xsz + i] = C(0.25)*(NW+NE+SW+SE);
				F1x[j*Parptr->xsz + i] = (NE+SE-NW-SW)/(C(4.0)*SQRT(C(3.0)));
				F1y[j*Parptr->xsz + i] = (NE+NW-SE-SW)/(C(4.0)*SQRT(C(3.0)));
			}
		}
	}
}

void dg2::downscale
(
	NUMERIC_TYPE* wd0,
	NUMERIC_TYPE* wd1x,
	NUMERIC_TYPE* wd1y,
	NUMERIC_TYPE no_data,
	Pars* Parptr,
	NUMERIC_TYPE* F
)
{
	NUMERIC_TYPE WD0, WD1x, WD1y;
	int ii, jj;
	//#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		for (int i = 0; i < Parptr->xsz; i++)
		{
			ii = 2 * i;
			jj = 4 * j;
			WD0 = wd0[j * Parptr->xsz + i];
			WD1x = wd1x[j * Parptr->xsz + i];
			WD1y = wd1y[j * Parptr->xsz + i];
			F[jj * Parptr->xsz + ii] = WD0 - C(0.5) * SQRT(C(3.0)) * WD1x + C(0.5) * SQRT(C(3.0)) * WD1y;
			F[jj * Parptr->xsz + ii + 1] = WD0 + C(0.5) * SQRT(C(3.0)) * WD1x + C(0.5) * SQRT(C(3.0)) * WD1y;
			F[jj * Parptr->xsz + ii + (2 * Parptr->xsz)] = WD0 - C(0.5) * SQRT(C(3.0)) * WD1x - C(0.5) * SQRT(C(3.0)) * WD1y;
			F[jj * Parptr->xsz + ii + (2 * Parptr->xsz) + 1] = WD0 + C(0.5) * SQRT(C(3.0)) * WD1x - C(0.5) * SQRT(C(3.0)) * WD1y;

		}
	}

	//#pragma omp parallel for private(ii, jj, WD0, WD1x, WD1y)
	//for (int j = 0; j < Parptr->ysz; j++)
	//{
	//	for (int i = 0; i < Parptr->xsz; i++)
	//	{
	//		ii = 4 * i;
	//		jj = 16 * j;
	//		WD0 = wd0[j * Parptr->xsz + i];
	//		WD1x = wd1x[j * Parptr->xsz + i];
	//		WD1y = wd1y[j * Parptr->xsz + i];
	//		F[jj * Parptr->xsz + ii]                          = WD0 - C(0.75) * SQRT(C(3.0)) * WD1x + C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + 1]                      = WD0 - C(0.25) * SQRT(C(3.0)) * WD1x + C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + 2]                      = WD0 + C(0.25) * SQRT(C(3.0)) * WD1x + C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + 3]                      = WD0 + C(0.75) * SQRT(C(3.0)) * WD1x + C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (4 * Parptr->xsz)]      = WD0 - C(0.75) * SQRT(C(3.0)) * WD1x + C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (4 * Parptr->xsz) + 1]  = WD0 - C(0.25) * SQRT(C(3.0)) * WD1x + C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (4 * Parptr->xsz) + 2]  = WD0 + C(0.25) * SQRT(C(3.0)) * WD1x + C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (4 * Parptr->xsz) + 3]  = WD0 + C(0.75) * SQRT(C(3.0)) * WD1x + C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (8 * Parptr->xsz)]      = WD0 - C(0.75) * SQRT(C(3.0)) * WD1x - C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (8 * Parptr->xsz) + 1]  = WD0 - C(0.25) * SQRT(C(3.0)) * WD1x - C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (8 * Parptr->xsz) + 2]  = WD0 + C(0.25) * SQRT(C(3.0)) * WD1x - C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (8 * Parptr->xsz) + 3]  = WD0 + C(0.75) * SQRT(C(3.0)) * WD1x - C(0.25) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (12 * Parptr->xsz)]     = WD0 - C(0.75) * SQRT(C(3.0)) * WD1x - C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (12 * Parptr->xsz) + 1] = WD0 - C(0.25) * SQRT(C(3.0)) * WD1x - C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (12 * Parptr->xsz) + 2] = WD0 + C(0.25) * SQRT(C(3.0)) * WD1x - C(0.75) * SQRT(C(3.0)) * WD1y;
	//		F[jj * Parptr->xsz + ii + (12 * Parptr->xsz) + 3] = WD0 + C(0.75) * SQRT(C(3.0)) * WD1x - C(0.75) * SQRT(C(3.0)) * WD1y;
	//	}
	//}
}

void dg2::initialise_h_from_eta
(
	NUMERIC_TYPE eta,
	Pars *Parptr,
	Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
			NUMERIC_TYPE Z1x = Arrptr->DEM1x[j*Parptr->xsz + i];
			NUMERIC_TYPE Z1y = Arrptr->DEM1y[j*Parptr->xsz + i];

			NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE& H1x = Arrptr->H1x[j*Parptr->xsz + i];
			NUMERIC_TYPE& H1y = Arrptr->H1y[j*Parptr->xsz + i];

			H = FMAX(C(0.0), eta - Z);
			H1x = -Z1x;
			H1y = -Z1y;
		}
	}
}

void dg2::read_dem_slopes
(
	Fnames *Fnameptr,
	Pars *Parptr,
	Arrays *Arrptr,
	const int verbose
)
{
	char dem1x[256];
	char dem1y[256];
	strcpy(dem1x, Fnameptr->demfilename);
	strcat(dem1x, "1x");
	strcpy(dem1y, Fnameptr->demfilename);
	strcat(dem1y, "1y");
	read_ascfile(dem1x, Parptr, Arrptr->DEM1x, "Loading DEM1x\n", verbose);
	read_ascfile(dem1y, Parptr, Arrptr->DEM1y, "Loading DEM1y\n", verbose);
}

void dg2::zero_dem_perimeter_slopes
(
    Pars *Parptr,
    Arrays *Arrptr
)
{
	for (int j=0; j<Parptr->ysz; j++)
	{
		// west
		{
			const int i=0;
			Arrptr->DEM1x[j*Parptr->xsz + i] = 0;
			Arrptr->DEM1y[j * Parptr->xsz + i] = 0;
		}

		// east
		{
			const int i=Parptr->xsz-1;
			Arrptr->DEM1x[j*Parptr->xsz + i] = 0;
			Arrptr->DEM1y[j * Parptr->xsz + i] = 0;
		}
	}

	for (int i=0; i<Parptr->xsz; i++)
	{
		// north
		{
			const int j=0;
			Arrptr->DEM1y[j*Parptr->xsz + i] = 0;
			Arrptr->DEM1x[j * Parptr->xsz + i] = 0;
		}

		// south
		{
			const int j=Parptr->ysz-1;
			Arrptr->DEM1y[j*Parptr->xsz + i] = 0;
			Arrptr->DEM1x[j * Parptr->xsz + i] = 0;
		}
	}
}

void dg2::read_h_slopes
(
	Fnames *Fnameptr,
	Pars *Parptr,
	Arrays *Arrptr,
	const int verbose
)
{
	char h1x[256];
	char h1y[256];
	strcpy(h1x, Fnameptr->startfilename);
	strcat(h1x, "1x");
	strcpy(h1y, Fnameptr->startfilename);
	strcat(h1y, "1y");
	read_ascfile(h1x, Parptr, Arrptr->H1x, "Loading startfile 1x\n", verbose);
	read_ascfile(h1y, Parptr, Arrptr->H1y, "Loading startfile 1y\n", verbose);
}

void dg2::read_discharge_slopes
(
	Fnames *Fnameptr,
	Pars *Parptr,
	Arrays *Arrptr,
	const int verbose
)
{
	char hu1x[256];
	char hu1y[256];
	char hv1x[256];
	char hv1y[256];
	strcpy(hu1x, Fnameptr->startfilename);
	strcat(hu1x, ".Qx1x");
	strcpy(hu1y, Fnameptr->startfilename);
	strcat(hu1y, ".Qx1y");
	strcpy(hv1x, Fnameptr->startfilename);
	strcat(hv1x, ".Qy1x");
	strcpy(hv1y, Fnameptr->startfilename);
	strcat(hv1y, ".Qy1y");
	read_ascfile(hu1x, Parptr, Arrptr->HU1x, "Loading startfile Qx1x\n", verbose);
	read_ascfile(hu1y, Parptr, Arrptr->HU1y, "Loading startfile Qx1y\n", verbose);
	read_ascfile(hv1x, Parptr, Arrptr->HV1x, "Loading startfile Qy1x\n", verbose);
	read_ascfile(hv1y, Parptr, Arrptr->HV1y, "Loading startfile Qy1y\n", verbose);
}

void dg2::write_dem
(
	const char *root,
	Pars *Parptr,
	Arrays *Arrptr
)
{
	States states;
	states.alt_ascheader = OFF;
	states.call_gzip = OFF;

	write_ascfile(root, -1, "", Arrptr->DEM, nullptr, 0,
			&states, Parptr, C(-1.0), "%.7");
	write_ascfile(root, -1, "1x", Arrptr->DEM1x, nullptr, 0,
			&states, Parptr, C(-1.0), "%.7");
	write_ascfile(root, -1, "1y", Arrptr->DEM1y, nullptr, 0,
			&states, Parptr, C(-1.0), "%.7");
}

void dg2::write_DScaled
(
	const char* root,
	Pars* Parptr,
	Arrays* Arrptr
)
{
	States states;
	states.alt_ascheader = OFF;
	states.call_gzip = OFF;

	write_ascfileDScaled(root, -1, "DScaled", Arrptr->H, nullptr, 0,
		&states, Parptr, C(-1.0), "%.7");
}

void dg2::write_startfile
(
	const char *root,
	Pars *Parptr,
	Arrays *Arrptr
)
{
	States states;
	states.alt_ascheader = OFF;
	states.call_gzip = OFF;

	write_ascfile(root, -1, "", Arrptr->H, nullptr, 0, &states, Parptr, C(-9999.0), "%.7");
	write_ascfile(root, -1, "1x", Arrptr->H1x, nullptr, 0, &states, Parptr, C(-9999.0), "%.7");
	write_ascfile(root, -1, "1y", Arrptr->H1y, nullptr, 0, &states, Parptr, C(-9999.0), "%.7");
}

void dg2::allocate_fields
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
	allocate_swe_fields(Parptr, Arrptr);
	allocate(Arrptr->boundary, Parptr);

	Arrptr->DEM1x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->DEM1y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->H1x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->H1y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HU1x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HU1y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV1x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV1y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);

	Arrptr->H_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->H1x_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->H1y_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HU_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HU1x_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HU1y_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV1x_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HV1y_int = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);

    if (Parptr->limit_slopes == ON)
    {
        Arrptr->ETA1x_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
        Arrptr->ETA1y_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
        Arrptr->HU1x_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
        Arrptr->HU1y_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
        Arrptr->HV1x_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
        Arrptr->HV1y_slopelim = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
    }
	
	Arrptr->HUstar_neg_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HUstar_pos_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HUstar_neg_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HUstar_pos_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HVstar_neg_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HVstar_pos_x = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HVstar_neg_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	Arrptr->HVstar_pos_y = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
}

void dg2::allocate
(
	BoundaryValues& boundary,
	Pars *Parptr
)
{
	boundary.DEM = memory_allocate_zero_numeric_legacy(2*Parptr->xsz + 2*Parptr->ysz);
	boundary.H = memory_allocate_zero_numeric_legacy(2*Parptr->xsz + 2*Parptr->ysz);
	boundary.HU = memory_allocate_zero_numeric_legacy(2*Parptr->xsz + 2*Parptr->ysz);
	boundary.HV = memory_allocate_zero_numeric_legacy(2*Parptr->xsz + 2*Parptr->ysz);
}

void dg2::deallocate_fields
(
    Pars *Parptr,
    Arrays *Arrptr
)
{
	deallocate_swe_fields(Arrptr);
	deallocate(Arrptr->boundary);

	memory_free_legacy(&Arrptr->DEM1x);
	memory_free_legacy(&Arrptr->DEM1y);
	memory_free_legacy(&Arrptr->H1x);
	memory_free_legacy(&Arrptr->H1y);
	memory_free_legacy(&Arrptr->HU1x);
	memory_free_legacy(&Arrptr->HU1y);
	memory_free_legacy(&Arrptr->HV1x);
	memory_free_legacy(&Arrptr->HV1y);

	memory_free_legacy(&Arrptr->H_int);
	memory_free_legacy(&Arrptr->H1x_int);
	memory_free_legacy(&Arrptr->H1y_int);
	memory_free_legacy(&Arrptr->HU_int);
	memory_free_legacy(&Arrptr->HU1x_int);
	memory_free_legacy(&Arrptr->HU1y_int);
	memory_free_legacy(&Arrptr->HV_int);
	memory_free_legacy(&Arrptr->HV1x_int);
	memory_free_legacy(&Arrptr->HV1y_int);

    if (Parptr->limit_slopes == ON)
    {
        memory_free_legacy(&Arrptr->ETA1x_slopelim);
        memory_free_legacy(&Arrptr->ETA1y_slopelim);
        memory_free_legacy(&Arrptr->HU1x_slopelim);
        memory_free_legacy(&Arrptr->HU1y_slopelim);
        memory_free_legacy(&Arrptr->HV1x_slopelim);
        memory_free_legacy(&Arrptr->HV1y_slopelim);
    }

	memory_free_legacy(&Arrptr->HUstar_neg_x);
	memory_free_legacy(&Arrptr->HUstar_pos_x);
	memory_free_legacy(&Arrptr->HUstar_neg_y);
	memory_free_legacy(&Arrptr->HUstar_pos_y);
	memory_free_legacy(&Arrptr->HVstar_neg_x);
	memory_free_legacy(&Arrptr->HVstar_pos_x);
	memory_free_legacy(&Arrptr->HVstar_neg_y);
	memory_free_legacy(&Arrptr->HVstar_pos_y);
}

void dg2::deallocate(BoundaryValues& boundary)
{
	memory_free_legacy(&boundary.DEM);
	memory_free_legacy(&boundary.H);
	memory_free_legacy(&boundary.HU);
	memory_free_legacy(&boundary.HV);
}
