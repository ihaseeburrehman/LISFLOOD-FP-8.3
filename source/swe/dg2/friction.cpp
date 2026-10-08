#include "friction.h"
#include "modifiedvars.h"

void dg2::apply_friction
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

            NUMERIC_TYPE n = (Arrptr->Manningsn == nullptr)
                ? Parptr->FPn : Arrptr->Manningsn[j*Parptr->xsz + i];

          	if (H <= Solverptr->DepthThresh) continue;

			friction(Solverptr, H, HU, HV, n);
            
                NUMERIC_TYPE H_x_lower = gauss_lower(Arrptr->H, Arrptr->H1x,
                        Parptr, i, j);
                NUMERIC_TYPE H_x_upper = gauss_upper(Arrptr->H, Arrptr->H1x,
                        Parptr, i, j);
                NUMERIC_TYPE HU_x_lower = gauss_lower(Arrptr->HU, Arrptr->HU1x,
                        Parptr, i, j);
                NUMERIC_TYPE HU_x_upper = gauss_upper(Arrptr->HU, Arrptr->HU1x,
                        Parptr, i, j);
                NUMERIC_TYPE HV_x_lower = gauss_lower(Arrptr->HV, Arrptr->HV1x,
                        Parptr, i, j);
                NUMERIC_TYPE HV_x_upper = gauss_upper(Arrptr->HV, Arrptr->HV1x,
                        Parptr, i, j);
                NUMERIC_TYPE H_y_lower = gauss_lower(Arrptr->H, Arrptr->H1y,
                        Parptr, i, j);
                NUMERIC_TYPE H_y_upper = gauss_upper(Arrptr->H, Arrptr->H1y,
                        Parptr, i, j);
                NUMERIC_TYPE HU_y_lower = gauss_lower(Arrptr->HU, Arrptr->HU1y,
                        Parptr, i, j);
                NUMERIC_TYPE HU_y_upper = gauss_upper(Arrptr->HU, Arrptr->HU1y,
                        Parptr, i, j);
                NUMERIC_TYPE HV_y_lower = gauss_lower(Arrptr->HV, Arrptr->HV1y,
                        Parptr, i, j);
                NUMERIC_TYPE HV_y_upper = gauss_upper(Arrptr->HV, Arrptr->HV1y,
                        Parptr, i, j);

                NUMERIC_TYPE& HU1x = Arrptr->HU1x[j*Parptr->xsz + i];
                NUMERIC_TYPE& HU1y = Arrptr->HU1y[j*Parptr->xsz + i];
                NUMERIC_TYPE& HV1x = Arrptr->HV1x[j*Parptr->xsz + i];
                NUMERIC_TYPE& HV1y = Arrptr->HV1y[j*Parptr->xsz + i];

                friction(Solverptr, H_x_lower, HU_x_lower, HV_x_lower, n);
                friction(Solverptr, H_x_upper, HU_x_upper, HV_x_upper, n);
                friction(Solverptr, H_y_lower, HU_y_lower, HV_y_lower, n);
                friction(Solverptr, H_y_upper, HU_y_upper, HV_y_upper, n);

                HU1x = C(0.5)*(HU_x_upper - HU_x_lower);
                HU1y = C(0.5)*(HU_y_upper - HU_y_lower);
                HV1x = C(0.5)*(HV_x_upper - HV_x_lower);
                HV1y = C(0.5)*(HV_y_upper - HV_y_lower);
		}
	}
}

void dg2::friction
(
	Solver *Solverptr,
	NUMERIC_TYPE H,
	NUMERIC_TYPE& HU,
	NUMERIC_TYPE& HV,
	NUMERIC_TYPE n
)
{
	if (H <= Solverptr->DepthThresh)
	{
		HU = C(0.0);
		HV = C(0.0);
		return;
	}

	NUMERIC_TYPE U = HU/H;
	NUMERIC_TYPE V = HV/H;

	if (FABS(U) <= Solverptr->SpeedThresh && FABS(V) <= Solverptr->SpeedThresh)
	{
		return;
	}

	NUMERIC_TYPE speed = SQRT(U*U+V*V);
	NUMERIC_TYPE Cf = Solverptr->g*n*n / pow(H, C(1.0)/C(3.0));
	NUMERIC_TYPE Sfx = -Cf*U*speed;
	NUMERIC_TYPE Sfy = -Cf*V*speed;
	NUMERIC_TYPE Dx = C(1.0) + Solverptr->Tstep*Cf/H * (C(2.0)*U*U + V*V)/speed;
	NUMERIC_TYPE Dy = C(1.0) + Solverptr->Tstep*Cf/H * (U*U + C(2.0)*V*V)/speed;
	HU += Solverptr->Tstep * Sfx / Dx;
	HV += Solverptr->Tstep * Sfy / Dy;
}
