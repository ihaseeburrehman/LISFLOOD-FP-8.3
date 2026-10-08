#include "dg2.h"
#include "dg2/dg2_output.h"
#include "dg2/fields.h"
#include "dg2/friction.h"
#include "dg2/modifiedvars.h"
#include "dg2/slope_limiter.h"
#include "stats.h"
#include "boundary.h"
#include "flux.h"
#include "hll.h"
#include "output.h"
#include <algorithm>

void dg2::solve
(
	Fnames *Fnameptr,
	Files *Fptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Stage *Stageptr,
	Arrays *Arrptr,
	const int verbose
)
{
	read_dem_slopes(Fnameptr, Parptr, Arrptr, verbose);
    zero_dem_perimeter_slopes(Parptr, Arrptr);

	if (Statesptr->startfile == ON)
	{
		read_h_slopes(Fnameptr, Parptr, Arrptr, verbose);
		if (Statesptr->startq2d == ON)
		{
			read_discharge_slopes(Fnameptr, Parptr, Arrptr, verbose);
		}
	}

//    zero_thin_depth_slopes(Parptr, Solverptr, Arrptr); // removed in new code

	if (Statesptr->adaptive_ts == ON)
	{
		Solverptr->Tstep = 1e-6;
	}
	else
	{
		Solverptr->Tstep = Solverptr->InitTstep;
	}
	Solverptr->MinTstep = Solverptr->InitTstep;
	Solverptr->vol1 = DomainVol(Statesptr, Parptr, nullptr, Arrptr, nullptr);
	initialise_Zstar(Parptr, Arrptr); // removed in new code
    DynamicRain<> rain(Fnameptr->dynamicrainfilename, verbose);

	time_t loop_start;
	time(&loop_start);

	while (Solverptr->t < Solverptr->Sim_Time)
	{	

		zero_flux_stats(BCptr);

        if (Parptr->limit_slopes == ON) Solverptr->maxH = maxH(Parptr, Arrptr);
        rain.update_H(Parptr, Solverptr, Arrptr);
		update_point_sources(Parptr, Solverptr, BCptr, Arrptr);
		update_boundary_values(Parptr, Solverptr, BCptr, Arrptr);

		if (Statesptr->adaptive_ts == ON)
		{
			Solverptr->Tstep = Tstep_from_cfl(Parptr, Solverptr, Arrptr);
		}

		if (verbose == ON) printf("t=%f\tdt=%f\n", Solverptr->t, Solverptr->Tstep);

		apply_friction(Parptr, Solverptr, Arrptr);
		rk_stage1(Parptr, Solverptr, BCptr, Arrptr);
		rk_stage2(Parptr, Solverptr, BCptr, Arrptr);

        update_max_field(Parptr, Arrptr);
		accumulate_boundary_flux_stats(Parptr, Solverptr, BCptr, Arrptr);

		if (Solverptr->t > C(0.0))
		{
			Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
		}
		Solverptr->t += Solverptr->Tstep;
		Solverptr->itCount += 1;

		if (Solverptr->t >= Parptr->MassTotal) {
			update_mass_stats(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
			write_mass_stats(Fptr, Parptr, Solverptr, BCptr);

            if (Statesptr->save_stages == ON)
            {
                write_depth_samples(Fptr, Parptr, Solverptr, Stageptr, Arrptr);
                if (Statesptr->voutput_stage == ON)
                {
                    write_speed_samples(Fptr, Parptr, Solverptr, Stageptr, Arrptr);
                }
            }
		}

		if (Solverptr->t >= Parptr->SaveTotal) {
			if (Statesptr->voutput == ON)
			{
				update_velocity(Parptr, Solverptr, Arrptr);
			}
			write_solution_slopes(Fnameptr, Statesptr, Parptr, Solverptr,
					Arrptr);
			write_solution(Fnameptr, Statesptr, Parptr, Solverptr, Arrptr);
		}
	}

	time_t loop_end;
	time(&loop_end);

	double seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);

    write_max_field(Fnameptr, Statesptr, Parptr, Arrptr);

	if (verbose == ON) printf("Finished.\n\n");
}

void dg2::update_point_sources
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	for (int ps_i = 0; ps_i < BCptr->numPS; ps_i++)
	{
		const int i = BCptr->xpi[ps_i];
		const int j = BCptr->ypi[ps_i];
		
		NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];

		switch (BCptr->PS_Ident[ps_i])
		{
		case HFIX2:
			{
				NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
				NUMERIC_TYPE H_new = min(C(0.0), BCptr->PS_Val[ps_i] - Z);
				NUMERIC_TYPE Q = (H_new - H)*Parptr->dA / Solverptr->Tstep;
				H = H_new;	

				accumulate_point_flux_stats(BCptr, Q);
			}
			break;
		case HVAR3:
			{
				NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
				NUMERIC_TYPE H_new = min(C(0.0), linear_interpolate(
						BCptr->PS_TimeSeries[ps_i], Solverptr->t) - Z);
				NUMERIC_TYPE Q = (H_new - H)*Parptr->dA / Solverptr->Tstep;
				H = H_new;	
				
				accumulate_point_flux_stats(BCptr, Q);
			}
			break;
		case QFIX4:
			{
				NUMERIC_TYPE Q = BCptr->PS_Val[ps_i];
				H += Q * Parptr->dx * Solverptr->Tstep / Parptr->dA;

				accumulate_point_flux_stats(BCptr, Q*Parptr->dx);
			}
			break;
		case QVAR5:
			{
				NUMERIC_TYPE Q = linear_interpolate(
						BCptr->PS_TimeSeries[ps_i], Solverptr->t);
				H += Q * Parptr->dx * Solverptr->Tstep / Parptr->dA;

				accumulate_point_flux_stats(BCptr, Q*Parptr->dx);
			}
			break;
		case NONE0:
		case FREE1:
		default:
			break;
		}
	}
}

void dg2::update_boundary_values
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	// west
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = limit_pos(Arrptr->H, Arrptr->H1x, Parptr, i, j);
		NUMERIC_TYPE HU_inside = limit_pos(Arrptr->HU, Arrptr->HU1x,
				Parptr, i, j);
		NUMERIC_TYPE HV_inside = limit_pos(Arrptr->HV, Arrptr->HV1x,
				Parptr, i, j);

		NUMERIC_TYPE H_const = Arrptr->H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = Arrptr->HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = Arrptr->HV[j*Parptr->xsz + i];

		const int bc_i = boundary_index_w(Parptr, i, j);
		NUMERIC_TYPE& H_outside = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE& HU_outside = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE& HV_outside = Arrptr->boundary.HV[bc_i];

		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, 1, Zstar,
				H_const, HU_const, HV_const,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);
	}
	
	// east
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz-1;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i+1];

		NUMERIC_TYPE H_inside = limit_neg(Arrptr->H, Arrptr->H1x, Parptr, i, j);
		NUMERIC_TYPE HU_inside = limit_neg(Arrptr->HU, Arrptr->HU1x,
				Parptr, i, j);
		NUMERIC_TYPE HV_inside = limit_neg(Arrptr->HV, Arrptr->HV1x,
				Parptr, i, j);

		NUMERIC_TYPE H_const = Arrptr->H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = Arrptr->HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = Arrptr->HV[j*Parptr->xsz + i];

		const int bc_i = boundary_index_e(Parptr, i, j);
		NUMERIC_TYPE& H_outside = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE& HU_outside = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE& HV_outside = Arrptr->boundary.HV[bc_i];

		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, -1, Zstar,
				H_const, HU_const, HV_const,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);
	}

	// north
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = limit_neg(Arrptr->H, Arrptr->H1y, Parptr, i, j);
		NUMERIC_TYPE HU_inside = limit_neg(Arrptr->HU, Arrptr->HU1y,
				Parptr, i, j);
		NUMERIC_TYPE HV_inside = limit_neg(Arrptr->HV, Arrptr->HV1y,
				Parptr, i, j);

		NUMERIC_TYPE H_const = Arrptr->H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = Arrptr->HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = Arrptr->HV[j*Parptr->xsz + i];

		const int bc_i = boundary_index_n(Parptr, i, j);
		NUMERIC_TYPE& H_outside = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE& HU_outside = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE& HV_outside = Arrptr->boundary.HV[bc_i];

		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, -1, Zstar,
				H_const, HV_const, HU_const,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);
	}

	// south
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz-1;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[(j+1)*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = limit_pos(Arrptr->H, Arrptr->H1y, Parptr, i, j);
		NUMERIC_TYPE HU_inside = limit_pos(Arrptr->HU, Arrptr->HU1y,
				Parptr, i, j);
		NUMERIC_TYPE HV_inside = limit_pos(Arrptr->HV, Arrptr->HV1y,
				Parptr, i, j);

		NUMERIC_TYPE H_const = Arrptr->H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = Arrptr->HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = Arrptr->HV[j*Parptr->xsz + i];

		const int bc_i = boundary_index_s(Parptr, i, j);
		NUMERIC_TYPE& H_outside = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE& HU_outside = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE& HV_outside = Arrptr->boundary.HV[bc_i];

		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, 1, Zstar,
				H_const, HV_const, HU_const,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);
	}
}

void dg2::rk_stage1
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	FlowCoefficients U;
	set_initial_coefficients(Arrptr, U);

	if (Parptr->limit_slopes == ON)
    {
        apply_slope_limiter(Parptr, Solverptr, Arrptr, U);
    }
    else
    {
        zero_perimeter_slopes(Parptr, Arrptr, U);
    }
	update_Hstar(Parptr, Solverptr, Arrptr, U);
	update_HUstar(Parptr, Solverptr, Arrptr, U);
	update_HVstar(Parptr, Solverptr, Arrptr, U);
	update_fluxes(Parptr, Solverptr, Arrptr, U);
	update_fluxes_on_boundaries(Parptr, Solverptr, BCptr, Arrptr, U);

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			Increment U_inc = {};
			L(Parptr, Solverptr, Arrptr, i, j, U, U_inc);

			update_intermediate_coefficients(
					Parptr, Solverptr, Arrptr, i, j, U_inc);

            zero_discharge(Parptr, Solverptr, i, j,
                    Arrptr->H_int, Arrptr->H1x_int, Arrptr->H1y_int,
                    Arrptr->HU_int, Arrptr->HU1x_int, Arrptr->HU1y_int,
                    Arrptr->HV_int, Arrptr->HV1x_int, Arrptr->HV1y_int);
		}
	}
}

void dg2::rk_stage2
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	FlowCoefficients U_int;
	set_intermediate_coefficients(Arrptr, U_int);
	
	if (Parptr->limit_slopes == ON)
    {
        apply_slope_limiter(Parptr, Solverptr, Arrptr, U_int);
    }
    else
    {
        zero_perimeter_slopes(Parptr, Arrptr, U_int);
    }
	update_Hstar(Parptr, Solverptr, Arrptr, U_int);
	update_HUstar(Parptr, Solverptr, Arrptr, U_int);
	update_HVstar(Parptr, Solverptr, Arrptr, U_int);
	update_fluxes(Parptr, Solverptr, Arrptr, U_int);
	update_fluxes_on_boundaries(Parptr, Solverptr, BCptr, Arrptr, U_int);

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			Increment U_inc = {};
			L(Parptr, Solverptr, Arrptr, i, j, U_int, U_inc);

			update_final_coefficients(
					Parptr, Solverptr, Arrptr, i, j, U_inc);

            zero_discharge(Parptr, Solverptr, i, j,
                    Arrptr->H, Arrptr->H1x, Arrptr->H1y,
                    Arrptr->HU, Arrptr->HU1x, Arrptr->HU1y,
                    Arrptr->HV, Arrptr->HV1x, Arrptr->HV1y);
		}
	}
}

void dg2::update_fluxes
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=1; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H_neg = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
			NUMERIC_TYPE HU_neg = Arrptr->HUstar_neg_x[j*Parptr->xsz + i-1];
			NUMERIC_TYPE HV_neg = Arrptr->HVstar_neg_x[j*Parptr->xsz + i-1];

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = Arrptr->HUstar_pos_x[j*Parptr->xsz + i];
			NUMERIC_TYPE HV_pos = Arrptr->HVstar_pos_x[j*Parptr->xsz + i];

			NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

			HLL_x(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
				FHx, FHUx, FHVx);
		}
	}

#pragma omp parallel for
	for (int j=1; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H_neg = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_neg = Arrptr->HUstar_neg_y[j*Parptr->xsz + i];
			NUMERIC_TYPE HV_neg = Arrptr->HVstar_neg_y[j*Parptr->xsz + i];

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = Arrptr->HUstar_pos_y[(j-1)*Parptr->xsz + i];
			NUMERIC_TYPE HV_pos = Arrptr->HVstar_pos_y[(j-1)*Parptr->xsz + i];

			NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

			HLL_y(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
					FHy, FHUy, FHVy);
		}
	}
}

void dg2::update_fluxes_on_boundaries
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr,
	FlowCoefficients const& U
)
{
	// west
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = Arrptr->HUstar_pos_x[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_inside = Arrptr->HVstar_pos_x[j*Parptr->xsz + i];

		NUMERIC_TYPE& H_const = U.H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = U.HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = U.HV[j*Parptr->xsz + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		const int bc_i = boundary_index_w(Parptr, i, j);
		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, 1, Zstar,
				H_const, HU_const, HV_const,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

		HLL_x(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHx, FHUx, FHVx);
	}

	// east
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HU_inside = Arrptr->HUstar_neg_x[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HV_inside = Arrptr->HVstar_neg_x[j*Parptr->xsz + i-1];

		NUMERIC_TYPE& H_const = U.H[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HU_const = U.HU[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HV_const = U.HV[j*Parptr->xsz + i-1];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		const int bc_i = boundary_index_e(Parptr, i-1, j);
		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, -1, Zstar,
				H_const, HU_const, HV_const,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

		HLL_x(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHx, FHUx, FHVx);
	}
	
	// north
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = Arrptr->HUstar_neg_y[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_inside = Arrptr->HVstar_neg_y[j*Parptr->xsz + i];

		NUMERIC_TYPE& H_const = U.H[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = U.HU[j*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = U.HV[j*Parptr->xsz + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		const int bc_i = boundary_index_n(Parptr, i, j);
		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, -1, Zstar,
				H_const, HV_const, HU_const,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		HLL_y(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHy, FHUy, FHVy);
	}

	// south
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = Arrptr->HUstar_pos_y[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HV_inside = Arrptr->HVstar_pos_y[(j-1)*Parptr->xsz + i];

		NUMERIC_TYPE& H_const = U.H[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HU_const = U.HU[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HV_const = U.HV[(j-1)*Parptr->xsz + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		const int bc_i = boundary_index_s(Parptr, i, j);
		set_boundary_values(Parptr, Solverptr, BCptr, bc_i, 1, Zstar,
				H_const, HV_const, HU_const,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		HLL_y(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHy, FHUy, FHVy);
	}
}

void dg2::set_boundary_values
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	const int bc_i,
	const int HU_sign,
	NUMERIC_TYPE Zstar,
	NUMERIC_TYPE& H_const,
	NUMERIC_TYPE HU_const,
	NUMERIC_TYPE HV_const,
	NUMERIC_TYPE& H_inside,
	NUMERIC_TYPE& HU_inside,
	NUMERIC_TYPE& HV_inside,
	NUMERIC_TYPE& H_outside,
	NUMERIC_TYPE& HU_outside,
	NUMERIC_TYPE& HV_outside
)
{
	H_outside = H_const;
	HU_outside = HU_const;
	HV_outside = HV_const;

	switch (BCptr->BC_Ident[bc_i])
	{
	case FREE1:

		if (H_const <= Solverptr->DepthThresh){
			H_const = C(0.0);
		    H_outside = C(0.0);
	    }
		break;
	case HFIX2:
		//H_inside = 
		H_outside = FMAX(C(0.0), BCptr->BC_Val[bc_i] - Zstar);
		break;
	case HVAR3:
		//H_inside = 
		H_outside = FMAX(C(0.0), linear_interpolate(
				BCptr->BC_TimeSeries[bc_i], Solverptr->t) - Zstar);
		break;
	case QFIX4:
		HU_inside = HU_outside = HU_sign * BCptr->BC_Val[bc_i];
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
		else
		{
			H_outside = H_inside;
		}
		break;
	case QVAR5:
		HU_inside = HU_outside = HU_sign * linear_interpolate(
				BCptr->BC_TimeSeries[bc_i], Solverptr->t);
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
		else
		{
			H_outside = H_inside;
		}
		break;
	case NONE0:
	default:
		H_outside = H_inside;
		HU_outside = -HU_inside;
		break;
	}
}

void dg2::update_intermediate_coefficients
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays* Arrptr,
	const int i,
	const int j,
	Increment const& U_inc
)
{
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H, Arrptr->H_int, U_inc.H);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H1x, Arrptr->H1x_int, U_inc.H1x);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H1y, Arrptr->H1y_int, U_inc.H1y);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU, Arrptr->HU_int, U_inc.HU);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU1x, Arrptr->HU1x_int, U_inc.HU1x);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU1y, Arrptr->HU1y_int, U_inc.HU1y);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV, Arrptr->HV_int, U_inc.HV);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV1x, Arrptr->HV1x_int, U_inc.HV1x);
		update_intermediate_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV1y, Arrptr->HV1y_int, U_inc.HV1y);
}

void dg2::update_intermediate_coefficient
(
	Pars *Parptr,
	Solver *Solverptr,
	const int i,
	const int j,
	NUMERIC_TYPE *field_n,
	NUMERIC_TYPE *field_int,
	NUMERIC_TYPE increment
)
{
	field_int[j*Parptr->xsz + i] = field_n[j*Parptr->xsz + i]
		+ Solverptr->Tstep * increment;
}

void dg2::update_final_coefficients
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays* Arrptr,
	const int i,
	const int j,
	Increment const& U_inc
)
{
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H, Arrptr->H_int, U_inc.H);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H1x, Arrptr->H1x_int, U_inc.H1x);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->H1y, Arrptr->H1y_int, U_inc.H1y);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU, Arrptr->HU_int, U_inc.HU);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU1x, Arrptr->HU1x_int, U_inc.HU1x);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HU1y, Arrptr->HU1y_int, U_inc.HU1y);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV, Arrptr->HV_int, U_inc.HV);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV1x, Arrptr->HV1x_int, U_inc.HV1x);
		update_final_coefficient(Parptr, Solverptr, i, j,
			Arrptr->HV1y, Arrptr->HV1y_int, U_inc.HV1y);
}

void dg2::update_final_coefficient
(
	Pars *Parptr,
	Solver *Solverptr,
	const int i,
	const int j,
	NUMERIC_TYPE *field_n,
	NUMERIC_TYPE *field_int,
	NUMERIC_TYPE increment
)
{
	field_n[j*Parptr->xsz + i] = C(0.5)*(field_n[j*Parptr->xsz + i]
		+ field_int[j*Parptr->xsz + i] + Solverptr->Tstep * increment);
}

void dg2::set_initial_coefficients
(
	Arrays *Arrptr,
	FlowCoefficients& U
)
{
	U.H = Arrptr->H;
	U.H1x = Arrptr->H1x;
	U.H1y = Arrptr->H1y;
	U.HU = Arrptr->HU;
	U.HU1x = Arrptr->HU1x;
	U.HU1y = Arrptr->HU1y;
	U.HV = Arrptr->HV;
	U.HV1x = Arrptr->HV1x;
	U.HV1y = Arrptr->HV1y;
}

void dg2::set_intermediate_coefficients
(
	Arrays *Arrptr,
	FlowCoefficients& U_int
)
{
	U_int.H = Arrptr->H_int;
	U_int.H1x = Arrptr->H1x_int;
	U_int.H1y = Arrptr->H1y_int;
	U_int.HU = Arrptr->HU_int;
	U_int.HU1x = Arrptr->HU1x_int;
	U_int.HU1y = Arrptr->HU1y_int;
	U_int.HV = Arrptr->HV_int;
	U_int.HV1x = Arrptr->HV1x_int;
	U_int.HV1y = Arrptr->HV1y_int;
}

void dg2::L
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j,
	FlowCoefficients const& U,
	Increment& U_inc
)
{
	L0(Parptr, Solverptr, Arrptr, i, j, U, U_inc);
    L1x(Parptr, Solverptr, Arrptr, i, j, U, U_inc);
    L1y(Parptr, Solverptr, Arrptr, i, j, U, U_inc);
}

void dg2::L0
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j,
	FlowCoefficients const& U,
	Increment& U_inc
)
{
	NUMERIC_TYPE H_w = Arrptr->FHx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE H_e = Arrptr->FHx[j*(Parptr->xsz+1) + i+1];
	NUMERIC_TYPE H_n = Arrptr->FHy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE H_s = Arrptr->FHy[(j+1)*(Parptr->xsz+1) + i];

	U_inc.H = -((H_e - H_w)/Parptr->dx + (H_n - H_s)/Parptr->dy);

	NUMERIC_TYPE HU_w = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HU_e = Arrptr->FHUx[j*(Parptr->xsz+1) + i+1];
	NUMERIC_TYPE HU_n = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HU_s = Arrptr->FHUy[(j+1)*(Parptr->xsz+1) + i];

	U_inc.HU = -((HU_e - HU_w)/Parptr->dx + (HU_n - HU_s)/Parptr->dy
			+ bed_source_0x(Parptr, Solverptr, Arrptr, U, i, j));

	NUMERIC_TYPE HV_w = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HV_e = Arrptr->FHVx[j*(Parptr->xsz+1) + i+1];
	NUMERIC_TYPE HV_n = Arrptr->FHVy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HV_s = Arrptr->FHVy[(j+1)*(Parptr->xsz+1) + i];

	U_inc.HV = -((HV_e - HV_w)/Parptr->dx + (HV_n - HV_s)/Parptr->dy
			+ bed_source_0y(Parptr, Solverptr, Arrptr, U, i, j));
}

void dg2::L1x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j,
	FlowCoefficients const& U,
	Increment& U_inc
)
{
	NUMERIC_TYPE H_w = Arrptr->FHx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE H_e = Arrptr->FHx[j*(Parptr->xsz+1) + i+1];
	NUMERIC_TYPE HU_w = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HU_e = Arrptr->FHUx[j*(Parptr->xsz+1) + i+1];
	NUMERIC_TYPE HV_w = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HV_e = Arrptr->FHVx[j*(Parptr->xsz+1) + i+1];

	NUMERIC_TYPE H = Hstar0x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE H1x = Hstar1x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HU = HUstar0x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HU1x = HUstar1x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HV = HVstar0x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HV1x = HVstar1x(Parptr, Arrptr, i, j);

	NUMERIC_TYPE H_g1 = gauss_lower(H, H1x);
	NUMERIC_TYPE HU_g1 = gauss_lower(HU, HU1x);
	NUMERIC_TYPE HV_g1 = gauss_lower(HV, HV1x);

	NUMERIC_TYPE H_g1_flux = C(0.0);
	NUMERIC_TYPE HU_g1_flux = C(0.0);
	NUMERIC_TYPE HV_g1_flux = C(0.0);

	physical_flux_x(
			Solverptr, H_g1, HU_g1, HV_g1, H_g1_flux, HU_g1_flux, HV_g1_flux);

	NUMERIC_TYPE H_g2 = gauss_upper(H, H1x);
	NUMERIC_TYPE HU_g2 = gauss_upper(HU, HU1x);
	NUMERIC_TYPE HV_g2 = gauss_upper(HV, HV1x);

	NUMERIC_TYPE H_g2_flux = C(0.0);
	NUMERIC_TYPE HU_g2_flux = C(0.0);
	NUMERIC_TYPE HV_g2_flux = C(0.0);

	physical_flux_x(
			Solverptr, H_g2, HU_g2, HV_g2, H_g2_flux, HU_g2_flux, HV_g2_flux);

	U_inc.H1x = -SQRT(C(3.0))/Parptr->dx*(H_w + H_e - H_g1_flux - H_g2_flux);
	U_inc.HU1x = -SQRT(C(3.0))/Parptr->dx*(HU_w + HU_e
			- HU_g1_flux - HU_g2_flux
			+ bed_source_1x(Parptr, Solverptr, Arrptr, U, i, j));
	U_inc.HV1x = -SQRT(C(3.0))/Parptr->dx*(HV_w + HV_e
			- HV_g1_flux - HV_g2_flux);
}

void dg2::L1y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j,
	FlowCoefficients const& U,
	Increment& U_inc
)
{
	NUMERIC_TYPE H_n = Arrptr->FHy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE H_s = Arrptr->FHy[(j+1)*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HU_n = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HU_s = Arrptr->FHUy[(j+1)*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HV_n = Arrptr->FHVy[j*(Parptr->xsz+1) + i];
	NUMERIC_TYPE HV_s = Arrptr->FHVy[(j+1)*(Parptr->xsz+1) + i];

	NUMERIC_TYPE H = Hstar0y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE H1y = Hstar1y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HU = HUstar0y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HU1y = HUstar1y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HV = HVstar0y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE HV1y = HVstar1y(Parptr, Arrptr, i, j);

	NUMERIC_TYPE H_g1 = gauss_lower(H, H1y);
	NUMERIC_TYPE HU_g1 = gauss_lower(HU, HU1y);
	NUMERIC_TYPE HV_g1 = gauss_lower(HV, HV1y);

	NUMERIC_TYPE H_g1_flux = C(0.0);
	NUMERIC_TYPE HU_g1_flux = C(0.0);
	NUMERIC_TYPE HV_g1_flux = C(0.0);

	physical_flux_y(
			Solverptr, H_g1, HU_g1, HV_g1, H_g1_flux, HU_g1_flux, HV_g1_flux);

	NUMERIC_TYPE H_g2 = gauss_upper(H, H1y);
	NUMERIC_TYPE HU_g2 = gauss_upper(HU, HU1y);
	NUMERIC_TYPE HV_g2 = gauss_upper(HV, HV1y);

	NUMERIC_TYPE H_g2_flux = C(0.0);
	NUMERIC_TYPE HU_g2_flux = C(0.0);
	NUMERIC_TYPE HV_g2_flux = C(0.0);

	physical_flux_y(
			Solverptr, H_g2, HU_g2, HV_g2, H_g2_flux, HU_g2_flux, HV_g2_flux);

	U_inc.H1y = -SQRT(C(3.0))/Parptr->dy*(H_s + H_n - H_g1_flux - H_g2_flux);
	U_inc.HU1y = -SQRT(C(3.0))/Parptr->dy*(HU_s + HU_n
			- HU_g1_flux - HU_g2_flux);
	U_inc.HV1y = -SQRT(C(3.0))/Parptr->dy*(HV_s + HV_n
			- HV_g1_flux - HV_g2_flux
			+ bed_source_1y(Parptr, Solverptr, Arrptr, U, i, j));
}

NUMERIC_TYPE dg2::bed_source_0x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar = Hstar0x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger = Zdagger1x(Parptr, Arrptr, U, i, j);

	return C(2.0)*SQRT(C(3.0)) * Solverptr->g * Hstar * Zdagger / Parptr->dx;
}

NUMERIC_TYPE dg2::bed_source_0y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar = Hstar0y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger = Zdagger1y(Parptr, Arrptr, U, i, j);

	return C(2.0)*SQRT(C(3.0)) * Solverptr->g * Hstar * Zdagger / Parptr->dy;
}

NUMERIC_TYPE dg2::bed_source_1x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar = Hstar1x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger = Zdagger1x(Parptr, Arrptr, U, i, j);

	return C(2.0) * Solverptr->g * Hstar * Zdagger;
}

NUMERIC_TYPE dg2::bed_source_1y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar = Hstar1y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger = Zdagger1y(Parptr, Arrptr, U, i, j);

	return C(2.0) * Solverptr->g * Hstar * Zdagger;
}

NUMERIC_TYPE dg2::maxH
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
	NUMERIC_TYPE H = C(0.0);

#ifndef _MSC_VER
#pragma omp parallel for reduction(max:H)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			H = FMAX(H, Arrptr->H[j*Parptr->xsz + i]);
		}
	}

	return H;
}

void dg2::zero_thin_depth_slopes
(
    Pars *Parptr,
    Solver *Solverptr,
    Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
            if (thin_depth(Solverptr, Arrptr->H[j*Parptr->xsz + i]))
            {
				Arrptr->H1x[j*Parptr->xsz + i] = C(0.0);
				Arrptr->H1y[j*Parptr->xsz + i] = C(0.0);
				Arrptr->HU1x[j*Parptr->xsz + i] = C(0.0);
				Arrptr->HU1y[j*Parptr->xsz + i] = C(0.0);
				Arrptr->HV1x[j*Parptr->xsz + i] = C(0.0);
				Arrptr->HV1y[j*Parptr->xsz + i] = C(0.0);
            }
		}
	}
}

bool dg2::thin_depth
(
    Solver *Solverptr,
    NUMERIC_TYPE H
)
{
//    return H > Solverptr->DepthThresh && H <= Solverptr->DG2DepthThresh;
	return H <= C(10.0) * Solverptr->DepthThresh; // thin layer depth is 10*thresh
}

NUMERIC_TYPE dg2::Tstep_from_cfl
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
	NUMERIC_TYPE dt = Solverptr->InitTstep;

#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
            NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];
            NUMERIC_TYPE HU = Arrptr->HU[j*Parptr->xsz + i];
            NUMERIC_TYPE HV = Arrptr->HV[j*Parptr->xsz + i];

            dt = min_dt(dt, H, HU, HV, Parptr, Solverptr);
		}
	}

// west
#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		const int bc_i = boundary_index_w(Parptr, i, j);
		NUMERIC_TYPE H = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE HU = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE HV = Arrptr->boundary.HV[bc_i];

		dt = min_dt(dt, H, HU, HV, Parptr, Solverptr);
	}

// east
#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz - 1;
		const int bc_i = boundary_index_e(Parptr, i, j);
		NUMERIC_TYPE H = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE HU = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE HV = Arrptr->boundary.HV[bc_i];

		dt = min_dt(dt, H, HU, HV, Parptr, Solverptr);
	}

// north
#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		const int bc_i = boundary_index_n(Parptr, i, j);
		NUMERIC_TYPE H = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE HU = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE HV = Arrptr->boundary.HV[bc_i];

		dt = min_dt(dt, H, HU, HV, Parptr, Solverptr);
	}

// south
#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz - 1;
		const int bc_i = boundary_index_s(Parptr, i, j);
		NUMERIC_TYPE H = Arrptr->boundary.H[bc_i];
		NUMERIC_TYPE HU = Arrptr->boundary.HU[bc_i];
		NUMERIC_TYPE HV = Arrptr->boundary.HV[bc_i];

		dt = min_dt(dt, H, HU, HV, Parptr, Solverptr);
	}

	return dt;
}
	
NUMERIC_TYPE dg2::min_dt
(
	NUMERIC_TYPE dt,
	NUMERIC_TYPE H,
	NUMERIC_TYPE HU,
	NUMERIC_TYPE HV,
    Pars *Parptr,
	Solver *Solverptr
)
{
    //if (thin_depth(Solverptr, H) && Solverptr->DG2ThinDepthTstep > C(0.0))
    //{
    //    return (std::min)({dt, Solverptr->DG2ThinDepthTstep});
    //}
    if (H > Solverptr->DepthThresh)
	{
		NUMERIC_TYPE U = HU/H;
		NUMERIC_TYPE V = HV/H;
		NUMERIC_TYPE dt_local_x = Solverptr->cfl*Parptr->dx/(FABS(U)+SQRT(Solverptr->g*H));
		NUMERIC_TYPE dt_local_y = Solverptr->cfl*Parptr->dy/(FABS(V)+SQRT(Solverptr->g*H));
		return (std::min)({dt, dt_local_x, dt_local_y});
	}
	else
	{
		return dt;
	}
}

void dg2::zero_discharge
(
    Pars* Parptr,
    Solver* Solverptr,
    int i,
    int j,
    NUMERIC_TYPE* H,
    NUMERIC_TYPE* H1x,
    NUMERIC_TYPE* H1y,
    NUMERIC_TYPE* HU,
    NUMERIC_TYPE* HU1x,
    NUMERIC_TYPE* HU1y,
    NUMERIC_TYPE* HV,
    NUMERIC_TYPE* HV1x,
    NUMERIC_TYPE* HV1y
)
{
    if (H[j*Parptr->xsz + i] < Solverptr->DepthThresh)
    {
        HU[j*Parptr->xsz + i] = C(0.0);
        HU1x[j*Parptr->xsz + i] = C(0.0);
        HU1y[j*Parptr->xsz + i] = C(0.0);
        HV[j*Parptr->xsz + i] = C(0.0);
        HV1x[j*Parptr->xsz + i] = C(0.0);
        HV1y[j*Parptr->xsz + i] = C(0.0);
    }
    //else if (thin_depth(Solverptr, H[j*Parptr->xsz + i]))
    //{
    //    H1x[j*Parptr->xsz + i] = C(0.0);
    //    H1y[j*Parptr->xsz + i] = C(0.0);
    //    HU1x[j*Parptr->xsz + i] = C(0.0);
    //    HU1y[j*Parptr->xsz + i] = C(0.0);
    //    HV1x[j*Parptr->xsz + i] = C(0.0);
    //    HV1y[j*Parptr->xsz + i] = C(0.0);
    //}
}
