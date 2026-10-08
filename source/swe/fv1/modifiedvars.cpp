#include "modifiedvars.h"

void fv1::initialise_Zstar
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
			NUMERIC_TYPE Z_neg = Arrptr->DEM[j*Parptr->xsz + i-1];
			NUMERIC_TYPE Z_pos = Arrptr->DEM[j*Parptr->xsz + i];
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
			NUMERIC_TYPE Z_neg = Arrptr->DEM[j*Parptr->xsz + i];
			NUMERIC_TYPE Z_pos = Arrptr->DEM[(j-1)*Parptr->xsz + i];
			NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];

			Zstar_y = getmax(Z_neg, Z_pos);
		}
	}

	// west boundary
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = 0;
		NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
		NUMERIC_TYPE& Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];
		Zstar_x = Z;
	}

	// east boundary
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		const int i = Parptr->xsz;
		NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i-1];
		NUMERIC_TYPE& Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];
		Zstar_x = Z;
	}

	// north boundary
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = 0;
		NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
		NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];
		Zstar_y = Z;
	}

	// south boundary
#pragma omp parallel for
	for (int i=0; i<Parptr->xsz; i++)
	{
		const int j = Parptr->ysz;
		NUMERIC_TYPE Z = Arrptr->DEM[(j-1)*Parptr->xsz + i];
		NUMERIC_TYPE& Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];
		Zstar_y = Z;
	}
}

void fv1::update_Hstar
(
	Pars *Parptr,
	Arrays *Arrptr
)
{
#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE ETA = eta(Parptr, Arrptr, i, j);
			NUMERIC_TYPE Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i+1];
			NUMERIC_TYPE& Hstar_neg_x = Arrptr->Hstar_neg_x[j*Parptr->xsz + i];

			Hstar_neg_x = getmax(C(0.0), ETA - Zstar_x);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE ETA = eta(Parptr, Arrptr, i, j);
			NUMERIC_TYPE Zstar_x = Arrptr->Zstar_x[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& Hstar_pos_x = Arrptr->Hstar_pos_x[j*Parptr->xsz + i];

			Hstar_pos_x = getmax(C(0.0), ETA - Zstar_x);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE ETA = eta(Parptr, Arrptr, i, j);
			NUMERIC_TYPE Zstar_y = Arrptr->Zstar_y[j*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& Hstar_neg_y = Arrptr->Hstar_neg_y[j*Parptr->xsz + i];

			Hstar_neg_y = getmax(C(0.0), ETA - Zstar_y);
		}
	}

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
			NUMERIC_TYPE ETA = eta(Parptr, Arrptr, i, j);
			NUMERIC_TYPE Zstar_y = Arrptr->Zstar_y[(j+1)*(Parptr->xsz+1) + i];
			NUMERIC_TYPE& Hstar_pos_y = Arrptr->Hstar_pos_y[j*Parptr->xsz + i];

			Hstar_pos_y = getmax(C(0.0), ETA - Zstar_y);
		}
	}
}

NUMERIC_TYPE fv1::HUstar_neg_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_neg_x, Arrptr->HU, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HUstar_pos_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_pos_x, Arrptr->HU, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HUstar_neg_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_neg_y, Arrptr->HU, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HUstar_pos_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_pos_y, Arrptr->HU, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HVstar_neg_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_neg_x, Arrptr->HV, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HVstar_pos_x
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_pos_x, Arrptr->HV, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HVstar_neg_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_neg_y, Arrptr->HV, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::HVstar_pos_y
(
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return discharge_star(
			Arrptr->Hstar_pos_y, Arrptr->HV, Parptr, Solverptr, Arrptr, i, j);
}

NUMERIC_TYPE fv1::discharge_star
(
	NUMERIC_TYPE *Hstar,
	NUMERIC_TYPE *discharge_component,
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE U = speed(
			discharge_component, Parptr, Solverptr, Arrptr, i, j);
	NUMERIC_TYPE Hstar_value = Hstar[j*Parptr->xsz + i];

	return Hstar_value * U;
}

NUMERIC_TYPE fv1::Zdagger_neg_x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return Zdagger(Parptr, Arrptr, Arrptr->Zstar_x, i+1, j, i, j);
}

NUMERIC_TYPE fv1::Zdagger_pos_x
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return Zdagger(Parptr, Arrptr, Arrptr->Zstar_x, i, j, i, j);
}

NUMERIC_TYPE fv1::Zdagger_neg_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return Zdagger(Parptr, Arrptr, Arrptr->Zstar_y, i, j, i, j);
}

NUMERIC_TYPE fv1::Zdagger_pos_y
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	return Zdagger(Parptr, Arrptr, Arrptr->Zstar_y, i, j+1, i, j);
}

NUMERIC_TYPE fv1::Zdagger
(
	Pars *Parptr,
	Arrays *Arrptr,
	NUMERIC_TYPE *Zstar,
	const int Zstar_i,
	const int Zstar_j,
	const int ETA_i,
	const int ETA_j
)
{
	NUMERIC_TYPE Zstar_value = Zstar[Zstar_j*(Parptr->xsz+1) + Zstar_i];
	NUMERIC_TYPE ETA = eta(Parptr, Arrptr, ETA_i, ETA_j);

	return Zstar_value - getmax(C(0.0), -(ETA - Zstar_value));
}

NUMERIC_TYPE fv1::eta
(
	Pars *Parptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE Z = Arrptr->DEM[j*Parptr->xsz + i];
	NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];

	return Z + H;
}

NUMERIC_TYPE fv1::speed
(
	NUMERIC_TYPE *discharge_component,
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr,
	const int i,
	const int j
)
{
	NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz + i];

	if (H > Solverptr->DepthThresh)
	{
		NUMERIC_TYPE discharge = discharge_component[j*Parptr->xsz + i];
		return discharge / H;
	}
	else
	{
		return C(0.0);
	}
}
