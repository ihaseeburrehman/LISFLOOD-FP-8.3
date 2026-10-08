#include "stats.h"

void update_mass_stats
(
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	Solverptr->FArea = flood_area(Parptr, Arrptr);
	Solverptr->vol2 = DomainVol(Statesptr, Parptr, nullptr, Arrptr, nullptr);
	Solverptr->Verror = BCptr->VolInMT - BCptr->VolOutMT -
		(Solverptr->vol2 - Solverptr->vol1);
	Solverptr->Qerror = Solverptr->Verror / Parptr->MassInt;

	BCptr->VolInMT = C(0.0);
	BCptr->VolOutMT = C(0.0);
	Solverptr->vol1 = Solverptr->vol2;
	Parptr->MassTotal += Parptr->MassInt;
}

NUMERIC_TYPE flood_area
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
	NUMERIC_TYPE FloodArea = C(0.0);
#pragma omp parallel for reduction(+:FloodArea)
	for (int j=0; j<Parptr->ysz; j++)
	{
		for (int i=0; i<Parptr->xsz; i++)
		{
			if (Arrptr->H[j*Parptr->xsz + i] > C(0.01)) FloodArea += Parptr->dA;
		}
	}
	return FloodArea;
}

void zero_flux_stats(BoundCs *BCptr)
{
	BCptr->Qin = C(0.0);
	BCptr->Qout = C(0.0);
}

void accumulate_point_flux_stats
(
	BoundCs *BCptr,
	NUMERIC_TYPE Q
)
{
	// no need to accumulate BCptr->Qpoint_pos or Qpoint_neg
	// just accumulate Qin and Qout directly
	if (Q > 0)
	{
		BCptr->Qin += Q;
	}
	else
	{
		BCptr->Qout -= Q;
	}
}

void accumulate_boundary_flux_stats
(
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr,
	Arrays *Arrptr
)
{
	NUMERIC_TYPE Qin = C(0.0);
	NUMERIC_TYPE Qout = C(0.0);

	// north
#pragma omp parallel for reduction(+:Qin, Qout)
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		if (FHy	< 0) Qin += -FHy*Parptr->dx; else Qout += FHy*Parptr->dx;
	}

	// south
#pragma omp parallel for reduction(+:Qin, Qout)
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE FHy = Arrptr->FHy[j*(Parptr->xsz+1) + i];
		if (FHy	> 0) Qin += FHy*Parptr->dx; else Qout += -FHy*Parptr->dx;
	}

	// west
#pragma omp parallel for reduction(+:Qin, Qout)
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		if (FHx	> 0) Qin += FHx*Parptr->dy; else Qout += -FHx*Parptr->dy;
	}

	// east
#pragma omp parallel for reduction(+:Qin, Qout)
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE FHx = Arrptr->FHx[j*(Parptr->xsz+1) + i];
		if (FHx	< 0) Qin += -FHx*Parptr->dy; else Qout += FHx*Parptr->dy;
	}

	BCptr->Qin += Qin;
	BCptr->Qout += Qout;
	BCptr->VolInMT += BCptr->Qin*Solverptr->Tstep;
	BCptr->VolOutMT += BCptr->Qout*Solverptr->Tstep;
}

void update_velocity
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
			NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];
			NUMERIC_TYPE HU = Arrptr->HU[j*Parptr->xsz + i];
			NUMERIC_TYPE HV = Arrptr->HV[j*Parptr->xsz + i];
			NUMERIC_TYPE& Vx = Arrptr->Vx[j*Parptr->xsz + i];
			NUMERIC_TYPE& Vy = Arrptr->Vy[j*Parptr->xsz + i];

			if (H > Solverptr->DepthThresh)
			{
				Vx = HU / H;
				Vy = HV / H;
			}
			else
			{
				Vx = C(0.0);
				Vy = C(0.0);
			}
		}
	}
}

void update_max_field
(
    Pars* Parptr,
    Arrays* Arrptr
)
{
#ifndef _MSC_VER
#pragma omp parallel for
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
            NUMERIC_TYPE& maxH = Arrptr->maxH[j*Parptr->xsz + i];
            maxH = FMAX(maxH, Arrptr->H[j*Parptr->xsz + i]);
        }
    }
}

