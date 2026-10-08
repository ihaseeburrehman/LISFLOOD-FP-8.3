#include "fv1.h"
#include "fv1/modifiedvars.h"
#include "boundary.h"
#include "hll.h"
#include "output.h"
#include "stats.h"
#include <algorithm>
#include <stdio.h>

void fv1::solve
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
	initialise_Zstar(Parptr, Arrptr);
    DynamicRain<> rain(Fnameptr->dynamicrainfilename, verbose);

	time_t loop_start;
	time(&loop_start);

	while (Solverptr->t < Solverptr->Sim_Time)
	{
		zero_flux_stats(BCptr);

        rain.update_H(Parptr, Solverptr, Arrptr);
		update_point_sources(Parptr, Solverptr, BCptr, Arrptr);
		apply_friction(Parptr, Solverptr, Arrptr);
		update_Hstar(Parptr, Arrptr);
		update_fluxes(Parptr, Solverptr, Arrptr);
		update_fluxes_on_boundaries(Parptr, Solverptr, BCptr, Arrptr);
		update_flow_variables(Parptr, Solverptr, Arrptr);
        fv1::drain_nodata_water(Parptr, Solverptr, BCptr, Arrptr);

		accumulate_boundary_flux_stats(Parptr, Solverptr, BCptr, Arrptr);

		if (Solverptr->t > C(0.0))
		{
			Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
		}
		Solverptr->t += Solverptr->Tstep;
		Solverptr->itCount += 1;

		if (Statesptr->adaptive_ts == ON)
		{
			Solverptr->Tstep = Tstep_from_cfl(Parptr, Solverptr, Arrptr);
		}

        update_max_field(Parptr, Arrptr);

		if (verbose == ON) printf("t=%f\tdt=%f\n", Solverptr->t, Solverptr->Tstep);

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

void fv1::apply_friction
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
	if (Arrptr->Manningsn == nullptr && Parptr->FPn <= C(0.0)) return;

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE& HU = Arrptr->HU[j*Parptr->xsz + i];
			NUMERIC_TYPE& HV = Arrptr->HV[j*Parptr->xsz + i];

			if (H <= Solverptr->DepthThresh) {
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

			NUMERIC_TYPE U = HU/H;
			NUMERIC_TYPE V = HV/H;
			if (FABS(U) <= Solverptr->SpeedThresh
					&& FABS(V) <= Solverptr->SpeedThresh)
			{
				HU = C(0.0);
				HV = C(0.0);
				continue;
			}

            NUMERIC_TYPE n = (Arrptr->Manningsn == nullptr)
                ? Parptr->FPn : Arrptr->Manningsn[j*Parptr->xsz + i];

			NUMERIC_TYPE Cf = Solverptr->g*n*n / pow(H, C(1.0)/C(3.0));
			NUMERIC_TYPE speed = SQRT(U*U+V*V);

			NUMERIC_TYPE Sf_x = -Cf*U*speed;
			NUMERIC_TYPE Sf_y = -Cf*V*speed;
			NUMERIC_TYPE D_x = C(1.0) + Solverptr->Tstep*Cf / H
				* (C(2.0)*U*U+V*V)/speed;
			NUMERIC_TYPE D_y = C(1.0) + Solverptr->Tstep*Cf / H
				* (U*U+C(2.0)*V*V)/speed;

			HU += Solverptr->Tstep*Sf_x/D_x;
			HV += Solverptr->Tstep*Sf_y/D_y;
		}
	}
}

void fv1::drain_nodata_water
(
    Pars* Parptr,
    Solver *Solverptr,
    BoundCs *BCptr,
    Arrays *Arrptr
)
{
    if (Parptr->drain_nodata == OFF) return;

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];

            if (FABS(Z - Parptr->nodata_elevation) < 1e-6)
            {
			    NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
                NUMERIC_TYPE Hold = H;
                H = C(0.0);

				accumulate_point_flux_stats(BCptr,
                        -Hold*Parptr->dA/Solverptr->Tstep);
            }
        }
    }
}

void fv1::update_point_sources
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
				NUMERIC_TYPE H_new = FMIN(C(0.0), BCptr->PS_Val[ps_i] - Z);
				NUMERIC_TYPE Q = (H_new - H)*Parptr->dA / Solverptr->Tstep;
				H = H_new;	

				accumulate_point_flux_stats(BCptr, Q);
			}
			break;
		case HVAR3:
			{
				NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
				NUMERIC_TYPE H_new = FMIN(C(0.0), linear_interpolate(
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

void fv1::update_fluxes
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=1; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE H_neg = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
			NUMERIC_TYPE HU_neg = HUstar_neg_x(
					Parptr, Solverptr, Arrptr, i-1, j);
			NUMERIC_TYPE HV_neg = HVstar_neg_x(
					Parptr, Solverptr, Arrptr, i-1, j);

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = HUstar_pos_x(Parptr, Solverptr, Arrptr, i, j);
			NUMERIC_TYPE HV_pos = HVstar_pos_x(Parptr, Solverptr, Arrptr, i, j);

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
			NUMERIC_TYPE HU_neg = HUstar_neg_y(Parptr, Solverptr, Arrptr, i, j);
			NUMERIC_TYPE HV_neg = HVstar_neg_y(Parptr, Solverptr, Arrptr, i, j);

			NUMERIC_TYPE H_pos = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
			NUMERIC_TYPE HU_pos = HUstar_pos_y(
					Parptr, Solverptr, Arrptr, i, j-1);
			NUMERIC_TYPE HV_pos = HVstar_pos_y(
					Parptr, Solverptr, Arrptr, i, j-1);

			NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

			HLL_y(Solverptr, H_neg, HU_neg, HV_neg, H_pos, HU_pos, HV_pos,
					FHy, FHUy, FHVy);
		}
	}
}

void fv1::update_fluxes_on_boundaries
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
		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_pos_x(Parptr, Solverptr, Arrptr, i, j);
		NUMERIC_TYPE HV_inside = HVstar_pos_x(Parptr, Solverptr, Arrptr, i, j);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
		
		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j,
				boundary_index_w(Parptr, i, j), 1,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		HLL_x(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHx, FHUx, FHVx);
	}

	// east
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_x[j*Parptr->xsz + i-1];
		NUMERIC_TYPE HU_inside = HUstar_neg_x(
				Parptr, Solverptr, Arrptr, i-1, j);
		NUMERIC_TYPE HV_inside = HVstar_neg_x(
				Parptr, Solverptr, Arrptr, i-1, j);

		NUMERIC_TYPE& FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUx = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVx = Arrptr->FHVx[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i-1, j,
				boundary_index_e(Parptr, i, j), -1,
				H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside);

		HLL_x(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHx, FHUx, FHVx);
	}
	
	// north
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_neg_y(Parptr, Solverptr, Arrptr, i, j);
		NUMERIC_TYPE HV_inside = HVstar_neg_y(Parptr, Solverptr, Arrptr, i, j);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j,
				boundary_index_n(Parptr, i, j), -1,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		HLL_y(Solverptr, H_inside, HU_inside, HV_inside,
				H_outside, HU_outside, HV_outside, FHy, FHUy, FHVy);
	}

	// south
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE H_inside = Arrptr->Hstar_pos_y[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE HU_inside = HUstar_pos_y(
				Parptr, Solverptr, Arrptr, i, j-1);
		NUMERIC_TYPE HV_inside = HVstar_pos_y(
				Parptr, Solverptr, Arrptr, i, j-1);

		NUMERIC_TYPE& FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHUy = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
		NUMERIC_TYPE& FHVy = Arrptr->FHVy[j*(Parptr->xsz+1) + i];

		NUMERIC_TYPE H_outside = C(0.0);
		NUMERIC_TYPE HU_outside = C(0.0);
		NUMERIC_TYPE HV_outside = C(0.0);

		set_boundary_values(Parptr, Solverptr, BCptr, Arrptr, i, j-1,
				boundary_index_s(Parptr, i, j), 1,
				H_inside, HV_inside, HU_inside,
				H_outside, HV_outside, HU_outside);

		HLL_y(Solverptr, H_outside, HU_outside, HV_outside,
				H_inside, HU_inside, HV_inside, FHy, FHUy, FHVy);
	}
}

void fv1::update_flow_variables
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
			NUMERIC_TYPE& H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE H_w = Arrptr->FHx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE H_e = Arrptr->FHx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE H_n = Arrptr->FHy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE H_s = Arrptr->FHy[(j+1)*(Parptr->xsz+1) + i];

			H = H - Solverptr->Tstep *
				((H_e - H_w)/Parptr->dx + (H_n - H_s)/Parptr->dy);

			NUMERIC_TYPE& HU = Arrptr->HU[j*Parptr->xsz + i];
			NUMERIC_TYPE HU_w = Arrptr->FHUx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HU_e = Arrptr->FHUx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE HU_n = Arrptr->FHUy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HU_s = Arrptr->FHUy[(j+1)*(Parptr->xsz+1) + i];

			HU = HU - Solverptr->Tstep *
				((HU_e - HU_w)/Parptr->dx + (HU_n - HU_s)/Parptr->dy
				 - bed_source_x(Parptr, Solverptr, Arrptr, i, j));

			NUMERIC_TYPE& HV = Arrptr->HV[j*Parptr->xsz + i];
			NUMERIC_TYPE HV_w = Arrptr->FHVx[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HV_e = Arrptr->FHVx[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE HV_n = Arrptr->FHVy[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE HV_s = Arrptr->FHVy[(j+1)*(Parptr->xsz+1) + i];

			HV = HV - Solverptr->Tstep *
				((HV_e - HV_w)/Parptr->dx + (HV_n - HV_s)/Parptr->dy
				 - bed_source_y(Parptr, Solverptr, Arrptr, i, j));
		}
	}
}

NUMERIC_TYPE fv1::bed_source_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg = Arrptr->Hstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_x(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_x(Parptr, Arrptr, i, j);

	return -Solverptr->g * C(0.5)*(Hstar_neg + Hstar_pos)
		* (Zdagger_neg - Zdagger_pos)/Parptr->dx;
}

NUMERIC_TYPE fv1::bed_source_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos = Arrptr->Hstar_pos_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_y(Parptr, Arrptr, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_y(Parptr, Arrptr, i, j);

	return -Solverptr->g * 0.5*(Hstar_neg + Hstar_pos)
		* (Zdagger_neg - Zdagger_pos)/Parptr->dy;
}

void fv1::set_boundary_values
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr,
	const int z_i,
	const int z_j,
	const int bc_i,
	const int HU_sign,
	NUMERIC_TYPE& H_inside,
	NUMERIC_TYPE& HU_inside,
	NUMERIC_TYPE& HV_inside,
	NUMERIC_TYPE& H_outside,
	NUMERIC_TYPE& HU_outside,
	NUMERIC_TYPE& HV_outside
)
{
	switch (BCptr->BC_Ident[bc_i])
	{
	case FREE1:
		H_outside = H_inside;
		HU_outside = HU_inside;
		HV_outside = HV_outside;
		break;
	case HFIX2:
		{
			NUMERIC_TYPE Z = Arrptr->DEM[z_j*Parptr->xsz + z_i];
			H_outside = FMAX(C(0.0), BCptr->BC_Val[bc_i] - Z);
			//H_inside = H_outside;
			HU_outside = HU_inside;
			HV_outside = HV_inside;
		}
		break;
	case HVAR3:
		{
			NUMERIC_TYPE Z = Arrptr->DEM[z_j*Parptr->xsz + z_i];
			H_outside = FMAX(C(0.0), linear_interpolate(
					BCptr->BC_TimeSeries[bc_i], Solverptr->t) - Z);
			//H_inside = H_outside;
			HU_outside = HU_inside;
			HV_outside = HV_inside;
		}
		break;
	case QFIX4:
		HU_outside = HU_sign * BCptr->BC_Val[bc_i];
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
        else
        {
		    H_outside = H_inside;
        }
		HU_inside = HU_outside;
		HV_outside = HV_inside;
		break;
	case QVAR5:
		HU_outside = HU_sign * linear_interpolate(
				BCptr->BC_TimeSeries[bc_i], Solverptr->t);
		if (FABS(HU_outside) > C(1e-12))
		{
			H_inside = H_outside = FMAX(H_inside, C(1.1)*Solverptr->DepthThresh);
		}
        else
        {
		    H_outside = H_inside;
        }
		HU_inside = HU_outside;
		HV_outside = HV_inside;
		break;
	case NONE0:
	default:
		H_outside = H_inside;
		HU_outside = -HU_inside;
		HV_outside = HV_outside;
		break;
	}
}

NUMERIC_TYPE fv1::Tstep_from_cfl
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
			if (H > Solverptr->DepthThresh)
			{
				NUMERIC_TYPE HU = Arrptr->HU[j*Parptr->xsz + i];
				NUMERIC_TYPE HV = Arrptr->HV[j*Parptr->xsz + i];
				NUMERIC_TYPE U = HU/H;
				NUMERIC_TYPE V = HV/H;
				
				NUMERIC_TYPE dt_x = Solverptr->cfl*Parptr->dx
					/ (FABS(U)+SQRT(Solverptr->g*H));
				NUMERIC_TYPE dt_y = Solverptr->cfl*Parptr->dy
					/ (FABS(V)+SQRT(Solverptr->g*H));

				dt = (std::min)({dt, dt_x, dt_y});
			}
		}
	}
	return dt;
}
