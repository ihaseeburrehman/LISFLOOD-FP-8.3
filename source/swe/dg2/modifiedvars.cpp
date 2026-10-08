#include "modifiedvars.h"

void dg2::initialise_Zstar
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
	// internal x
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=1; i<Parptr->xsz; i++)
		{
			
			NUMERIC_TYPE Z_neg = limit_neg(Arrptr->DEM, Arrptr->DEM1x,
					Parptr, i-1, j);

			NUMERIC_TYPE Z_pos = limit_pos(Arrptr->DEM, Arrptr->DEM1x,
					Parptr, i, j);

			NUMERIC_TYPE& Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];

			Zstar_x = getmax(Z_neg, Z_pos);
		}
	}

	// internal y
#pragma omp parallel for
	for (int j=1; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z_neg = limit_neg(Arrptr->DEM, Arrptr->DEM1y,
					Parptr, i, j);

			NUMERIC_TYPE Z_pos = limit_pos(Arrptr->DEM, Arrptr->DEM1y,
					Parptr, i, j-1);

			NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];

			Zstar_y = getmax(Z_neg, Z_pos);
		}
	}

	// west boundary
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE Z_pos = limit_pos(Arrptr->DEM, Arrptr->DEM1x, Parptr, i, j);
		NUMERIC_TYPE& Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];
		Zstar_x = Z_pos;
	}

	// east boundary
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE Z_neg = limit_neg(Arrptr->DEM, Arrptr->DEM1x, Parptr, i-1, j);
		NUMERIC_TYPE& Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];
		Zstar_x = Z_neg;
	}

	// north boundary
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE Z_pos = limit_pos(Arrptr->DEM, Arrptr->DEM1y, Parptr, i, j);
		NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];
		Zstar_y = Z_pos;
	}

	// south boundary
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE Z_neg = limit_neg(Arrptr->DEM, Arrptr->DEM1y, Parptr, i, j-1);
		NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];
		Zstar_y = Z_neg;
	}
}

void dg2::update_Hstar
(
	Pars *Parptr,
	Solver* Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
				NUMERIC_TYPE ETA_neg_x = eta_neg_x(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE Zstar_x = Arrptr->Zstar_x[j * (Parptr->xsz + 1) + i + 1];
				NUMERIC_TYPE& Hstar_neg_x = Arrptr->Hstar_neg_x[j * Parptr->xsz + i];

				Hstar_neg_x = getmax(C(0.0), ETA_neg_x - Zstar_x);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
				NUMERIC_TYPE ETA_pos_x = eta_pos_x(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE Zstar_x = Arrptr->Zstar_x[j * (Parptr->xsz + 1) + i];
				NUMERIC_TYPE& Hstar_pos_x = Arrptr->Hstar_pos_x[j * Parptr->xsz + i];

				Hstar_pos_x = getmax(C(0.0), ETA_pos_x - Zstar_x);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
				NUMERIC_TYPE ETA_neg_y = eta_neg_y(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE Zstar_y = Arrptr->Zstar_y[j * (Parptr->xsz + 1) + i];
				NUMERIC_TYPE& Hstar_neg_y = Arrptr->Hstar_neg_y[j * Parptr->xsz + i];

				Hstar_neg_y = getmax(C(0.0), ETA_neg_y - Zstar_y);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
				NUMERIC_TYPE ETA_pos_y = eta_pos_y(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE Zstar_y = Arrptr->Zstar_y[(j + 1) * (Parptr->xsz + 1) + i];
				NUMERIC_TYPE& Hstar_pos_y = Arrptr->Hstar_pos_y[j * Parptr->xsz + i];

				Hstar_pos_y = getmax(C(0.0), ETA_pos_y - Zstar_y);
		}
	}
}

void dg2::update_HUstar
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U
)
{
	update_discharge_star(Parptr, Solverptr, Arrptr, U.HU, U.HU1x, U.HU1y, U,
			Arrptr->HUstar_neg_x, Arrptr->HUstar_pos_x,
			Arrptr->HUstar_neg_y, Arrptr->HUstar_pos_y);
}

void dg2::update_HVstar
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	FlowCoefficients const& U
)
{
	update_discharge_star(Parptr, Solverptr, Arrptr, U.HV, U.HV1x, U.HV1y, U,
			Arrptr->HVstar_neg_x, Arrptr->HVstar_pos_x,
			Arrptr->HVstar_neg_y, Arrptr->HVstar_pos_y);
}

void dg2::update_discharge_star
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	NUMERIC_TYPE* discharge,
	NUMERIC_TYPE* discharge1x,
	NUMERIC_TYPE* discharge1y,
	FlowCoefficients const& U,
	NUMERIC_TYPE* discharge_star_neg_x,
	NUMERIC_TYPE* discharge_star_pos_x,
	NUMERIC_TYPE* discharge_star_neg_y,
	NUMERIC_TYPE* discharge_star_pos_y
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
				NUMERIC_TYPE H_neg_x = limit_neg(U.H, U.H1x, Parptr, i, j);
				NUMERIC_TYPE discharge_neg_x = limit_neg(discharge, discharge1x,
					Parptr, i, j);
				NUMERIC_TYPE Hstar_neg_x = Arrptr->Hstar_neg_x[j * Parptr->xsz + i];

				NUMERIC_TYPE H_pos_x = limit_pos(U.H, U.H1x, Parptr, i, j);
				NUMERIC_TYPE discharge_pos_x = limit_pos(discharge, discharge1x,
					Parptr, i, j);
				NUMERIC_TYPE Hstar_pos_x = Arrptr->Hstar_pos_x[j * Parptr->xsz + i];

				NUMERIC_TYPE H_neg_y = limit_neg(U.H, U.H1y, Parptr, i, j);
				NUMERIC_TYPE discharge_neg_y = limit_neg(discharge, discharge1y,
					Parptr, i, j);
				NUMERIC_TYPE Hstar_neg_y = Arrptr->Hstar_neg_y[j * Parptr->xsz + i];

				NUMERIC_TYPE H_pos_y = limit_pos(U.H, U.H1y, Parptr, i, j);
				NUMERIC_TYPE discharge_pos_y = limit_pos(discharge, discharge1y,
					Parptr, i, j);
				NUMERIC_TYPE Hstar_pos_y = Arrptr->Hstar_pos_y[j * Parptr->xsz + i];

				discharge_star_neg_x[j * Parptr->xsz + i] = discharge_star(Solverptr,
					H_neg_x, discharge_neg_x, Hstar_neg_x);
				discharge_star_pos_x[j * Parptr->xsz + i] = discharge_star(Solverptr,
					H_pos_x, discharge_pos_x, Hstar_pos_x);
				discharge_star_neg_y[j * Parptr->xsz + i] = discharge_star(Solverptr,
					H_neg_y, discharge_neg_y, Hstar_neg_y);
				discharge_star_pos_y[j * Parptr->xsz + i] = discharge_star(Solverptr,
					H_pos_y, discharge_pos_y, Hstar_pos_y);
		}
	}
}

NUMERIC_TYPE dg2::discharge_star
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE discharge,
	NUMERIC_TYPE Hstar
)
{
	if (H >= Solverptr->DepthThresh)
	{
		NUMERIC_TYPE speed = discharge / H;
		return Hstar * speed;
	}
	else
	{
		return C(0.0);
	}
}

NUMERIC_TYPE dg2::Hstar0x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg_x = Arrptr->Hstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos_x = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
	return C(0.5)*(Hstar_neg_x + Hstar_pos_x);
}

NUMERIC_TYPE dg2::Hstar0y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg_y = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos_y = Arrptr->Hstar_pos_y[j*Parptr->xsz + i];
	return C(0.5)*(Hstar_neg_y + Hstar_pos_y);
}

NUMERIC_TYPE dg2::Hstar1x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg_x = Arrptr->Hstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos_x = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];
	return (Hstar_neg_x - Hstar_pos_x)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::Hstar1y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Hstar_neg_y = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE Hstar_pos_y = Arrptr->Hstar_pos_y[j*Parptr->xsz + i];
	return (Hstar_neg_y - Hstar_pos_y)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::HUstar0x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HUstar_neg_x = Arrptr->HUstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE HUstar_pos_x = Arrptr->HUstar_pos_x[j*Parptr->xsz + i];
	return C(0.5)*(HUstar_neg_x + HUstar_pos_x);
}

NUMERIC_TYPE dg2::HUstar0y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HUstar_neg_y = Arrptr->HUstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE HUstar_pos_y = Arrptr->HUstar_pos_y[j*Parptr->xsz + i];
	return C(0.5)*(HUstar_neg_y + HUstar_pos_y);
}

NUMERIC_TYPE dg2::HUstar1x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HUstar_neg_x = Arrptr->HUstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE HUstar_pos_x = Arrptr->HUstar_pos_x[j*Parptr->xsz + i];
	return (HUstar_neg_x - HUstar_pos_x)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::HUstar1y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HUstar_neg_y = Arrptr->HUstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE HUstar_pos_y = Arrptr->HUstar_pos_y[j*Parptr->xsz + i];
	return (HUstar_neg_y - HUstar_pos_y)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::HVstar0x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HVstar_neg_x = Arrptr->HVstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE HVstar_pos_x = Arrptr->HVstar_pos_x[j*Parptr->xsz + i];
	return C(0.5)*(HVstar_neg_x + HVstar_pos_x);
}

NUMERIC_TYPE dg2::HVstar0y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HVstar_neg_y = Arrptr->HVstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE HVstar_pos_y = Arrptr->HVstar_pos_y[j*Parptr->xsz + i];
	return C(0.5)*(HVstar_neg_y + HVstar_pos_y);
}

NUMERIC_TYPE dg2::HVstar1x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HVstar_neg_x = Arrptr->HVstar_neg_x[j*Parptr->xsz + i];
	NUMERIC_TYPE HVstar_pos_x = Arrptr->HVstar_pos_x[j*Parptr->xsz + i];
	return (HVstar_neg_x - HVstar_pos_x)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::HVstar1y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE HVstar_neg_y = Arrptr->HVstar_neg_y[j*Parptr->xsz + i];
	NUMERIC_TYPE HVstar_pos_y = Arrptr->HVstar_pos_y[j*Parptr->xsz + i];
	return (HVstar_neg_y - HVstar_pos_y)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::Zdagger1x
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_x(Parptr, Arrptr, U, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_x(Parptr, Arrptr, U, i, j);

	return (Zdagger_neg - Zdagger_pos)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::Zdagger1y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Zdagger_neg = Zdagger_neg_y(Parptr, Arrptr, U, i, j);
	NUMERIC_TYPE Zdagger_pos = Zdagger_pos_y(Parptr, Arrptr, U, i, j);

	return (Zdagger_neg - Zdagger_pos)/(C(2.0)*SQRT(C(3.0)));
}

NUMERIC_TYPE dg2::Zdagger_neg_x
(
	Pars* Parptr,
	Arrays* Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{

		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j * (Parptr->xsz + 1) + i + 1];
		NUMERIC_TYPE ETA = eta_neg_x(Parptr, Arrptr, U, i, j);
		return Zstar - getmax(C(0.0), -(ETA - Zstar));


}

NUMERIC_TYPE dg2::Zdagger_pos_x
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{

		NUMERIC_TYPE Zstar = Arrptr->Zstar_x[j * (Parptr->xsz + 1) + i];
		NUMERIC_TYPE ETA = eta_pos_x(Parptr, Arrptr, U, i, j);

		return Zstar - getmax(C(0.0), -(ETA - Zstar));

}

NUMERIC_TYPE dg2::Zdagger_neg_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{

		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[j * (Parptr->xsz + 1) + i];
		NUMERIC_TYPE ETA = eta_neg_y(Parptr, Arrptr, U, i, j);

		return Zstar - getmax(C(0.0), -(ETA - Zstar));

}

NUMERIC_TYPE dg2::Zdagger_pos_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{

		NUMERIC_TYPE Zstar = Arrptr->Zstar_y[(j + 1) * (Parptr->xsz + 1) + i];
		NUMERIC_TYPE ETA = eta_pos_y(Parptr, Arrptr, U, i, j);

		return Zstar - getmax(C(0.0), -(ETA - Zstar));

}

NUMERIC_TYPE dg2::eta_neg_x
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE ETA = eta(Parptr, Arrptr, U, i, j);
    NUMERIC_TYPE DEM1x = Arrptr->DEM1x[j*Parptr->xsz + i];
    NUMERIC_TYPE H1x = U.H1x[j*Parptr->xsz + i];
	return limit_neg(ETA, DEM1x + H1x);
}

NUMERIC_TYPE dg2::eta_pos_x
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE ETA = eta(Parptr, Arrptr, U, i, j);
    NUMERIC_TYPE DEM1x = Arrptr->DEM1x[j*Parptr->xsz + i];
    NUMERIC_TYPE H1x = U.H1x[j*Parptr->xsz + i];
	return limit_pos(ETA, DEM1x + H1x);
}

NUMERIC_TYPE dg2::eta_neg_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE ETA = eta(Parptr, Arrptr, U, i, j);
    NUMERIC_TYPE DEM1y = Arrptr->DEM1y[j*Parptr->xsz + i];
    NUMERIC_TYPE H1y = U.H1y[j*Parptr->xsz + i];
	return limit_neg(ETA, DEM1y + H1y);
}

NUMERIC_TYPE dg2::eta_pos_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE ETA = eta(Parptr, Arrptr, U, i, j);
    NUMERIC_TYPE DEM1y = Arrptr->DEM1y[j*Parptr->xsz + i];
    NUMERIC_TYPE H1y = U.H1y[j*Parptr->xsz + i];
	return limit_pos(ETA, DEM1y + H1y);
}

NUMERIC_TYPE dg2::eta
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
	NUMERIC_TYPE H = U.H[j*Parptr->xsz + i];

	return Z + H;
}

NUMERIC_TYPE dg2::eta1x
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Z1x = Arrptr->DEM1x[j*Parptr->xsz + i];
	NUMERIC_TYPE H1x = U.H1x[j*Parptr->xsz + i];

	return Z1x + H1x;
}

NUMERIC_TYPE dg2::eta1y
(
	Pars *Parptr,
	Arrays *Arrptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Z1y = Arrptr->DEM1y[j*Parptr->xsz + i];
	NUMERIC_TYPE H1y = U.H1y[j*Parptr->xsz + i];

	return Z1y + H1y;
}

NUMERIC_TYPE dg2::limit_neg
(
	NUMERIC_TYPE *constant,
	NUMERIC_TYPE *slope,
	Pars *Parptr,
	const int i,
	const int j
)
{
	return limit_neg(constant[j*Parptr->xsz + i], slope[j*Parptr->xsz + i]);
}

NUMERIC_TYPE dg2::limit_pos
(
	NUMERIC_TYPE *constant,
	NUMERIC_TYPE *slope,
	Pars *Parptr,
	const int i,
	const int j
)
{
	return limit_pos(constant[j*Parptr->xsz + i], slope[j*Parptr->xsz + i]);
}

NUMERIC_TYPE dg2::gauss_lower
(
	NUMERIC_TYPE *constant,
	NUMERIC_TYPE *slope,
	Pars *Parptr,
	const int i,
	const int j
)
{
	return gauss_lower(constant[j*Parptr->xsz + i], slope[j*Parptr->xsz + i]);
}

NUMERIC_TYPE dg2::gauss_upper
(
	NUMERIC_TYPE *constant,
	NUMERIC_TYPE *slope,
	Pars *Parptr,
	const int i,
	const int j
)
{
	return gauss_upper(constant[j*Parptr->xsz + i], slope[j*Parptr->xsz + i]);
}

NUMERIC_TYPE dg2::limit_neg(NUMERIC_TYPE constant, NUMERIC_TYPE slope)
{
	return constant + SQRT(C(3.0))*slope;
}

NUMERIC_TYPE dg2::limit_pos(NUMERIC_TYPE constant, NUMERIC_TYPE slope)
{
	return constant - SQRT(C(3.0))*slope;
}

NUMERIC_TYPE dg2::gauss_lower(NUMERIC_TYPE constant, NUMERIC_TYPE slope)
{
	return constant - slope;
}

NUMERIC_TYPE dg2::gauss_upper(NUMERIC_TYPE constant, NUMERIC_TYPE slope)
{
	return constant + slope;
}
