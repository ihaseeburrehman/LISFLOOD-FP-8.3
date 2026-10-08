#include "slope_limiter.h"
#include "modifiedvars.h"
#include "../boundary.h"
#include <algorithm>

void dg2::zero_perimeter_slopes
(
	Pars* Parptr,
	Arrays* Arrptr,
	FlowCoefficients const& U
)
{
#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		// west
		{
			const int i = 0;
			U.H1x[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1x[j * Parptr->xsz + i];
			U.HU1x[j * Parptr->xsz + i] = C(0.0);
			U.HV1x[j * Parptr->xsz + i] = C(0.0);
			U.H1y[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1x[j * Parptr->xsz + i];
			U.HU1y[j * Parptr->xsz + i] = C(0.0);
			U.HV1y[j * Parptr->xsz + i] = C(0.0);
		}

		// east
		{
			const int i = Parptr->xsz - 1;
			U.H1x[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1x[j * Parptr->xsz + i];
			U.HU1x[j * Parptr->xsz + i] = C(0.0);
			U.HV1x[j * Parptr->xsz + i] = C(0.0);
			U.H1y[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1x[j * Parptr->xsz + i];
			U.HU1y[j * Parptr->xsz + i] = C(0.0);
			U.HV1y[j * Parptr->xsz + i] = C(0.0);
		}
	}

#pragma omp parallel for
	for (int i = 0; i < Parptr->xsz; i++)
	{
		// north
		{
			const int j = 0;
			U.H1y[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1y[j * Parptr->xsz + i];
			U.HU1y[j * Parptr->xsz + i] = C(0.0);
			U.HV1y[j * Parptr->xsz + i] = C(0.0);
			U.H1x[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1y[j * Parptr->xsz + i];
			U.HU1x[j * Parptr->xsz + i] = C(0.0);
			U.HV1x[j * Parptr->xsz + i] = C(0.0);
		}

		// south	
		{
			const int j = Parptr->ysz - 1;
			U.H1y[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1y[j * Parptr->xsz + i];
			U.HU1y[j * Parptr->xsz + i] = C(0.0);
			U.HV1y[j * Parptr->xsz + i] = C(0.0);
			U.H1x[j * Parptr->xsz + i] = C(0.0); //-Arrptr->DEM1y[j * Parptr->xsz + i];
			U.HU1x[j * Parptr->xsz + i] = C(0.0);
			U.HV1x[j * Parptr->xsz + i] = C(0.0);
		}
	}
}

void dg2::apply_slope_limiter
(
	Pars* Parptr,
	Solver* Solverptr,
	Arrays* Arrptr,
	FlowCoefficients const& U
)
{
	apply_slope_limiter_x(Parptr, Solverptr, Arrptr, U);
	apply_slope_limiter_y(Parptr, Solverptr, Arrptr, U);
}

void dg2::apply_slope_limiter_x
(
	Pars* Parptr,
	Solver* Solverptr,
	Arrays* Arrptr,
	FlowCoefficients const& U
)
{
#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		// west
		{
			const int i = 0;
			Arrptr->ETA1x_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HU1x_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HV1x_slopelim[j * Parptr->xsz + i] = C(0.0);
		}

		// east
		{
			const int i = Parptr->xsz - 1;
			Arrptr->ETA1x_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HU1x_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HV1x_slopelim[j * Parptr->xsz + i] = C(0.0);
		}
	}

#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		for (int i = 1; i < Parptr->xsz - 1; i++)
		{
			NUMERIC_TYPE minH = stencil_minH_x(Parptr, U, i, j);

			{
				NUMERIC_TYPE backward_const = eta(Parptr, Arrptr, U, i - 1, j);
				NUMERIC_TYPE backward_slope = eta1x(Parptr, Arrptr, U, i - 1, j);
				NUMERIC_TYPE local_const = eta(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE local_slope = eta1x(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE forward_const = eta(Parptr, Arrptr, U, i + 1, j);
				NUMERIC_TYPE forward_slope = eta1x(Parptr, Arrptr, U, i + 1, j);

				NUMERIC_TYPE& local_slopelim =
					Arrptr->ETA1x_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}

			{
				NUMERIC_TYPE backward_const = U.HU[j * Parptr->xsz + i - 1];
				NUMERIC_TYPE backward_slope = U.HU1x[j * Parptr->xsz + i - 1];
				NUMERIC_TYPE local_const = U.HU[j * Parptr->xsz + i];
				NUMERIC_TYPE local_slope = U.HU1x[j * Parptr->xsz + i];
				NUMERIC_TYPE forward_const = U.HU[j * Parptr->xsz + i + 1];
				NUMERIC_TYPE forward_slope = U.HU1x[j * Parptr->xsz + i + 1];

				NUMERIC_TYPE& local_slopelim =
					Arrptr->HU1x_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}

			{
				NUMERIC_TYPE backward_const = U.HV[j * Parptr->xsz + i - 1];
				NUMERIC_TYPE backward_slope = U.HV1x[j * Parptr->xsz + i - 1];
				NUMERIC_TYPE local_const = U.HV[j * Parptr->xsz + i];
				NUMERIC_TYPE local_slope = U.HV1x[j * Parptr->xsz + i];
				NUMERIC_TYPE forward_const = U.HV[j * Parptr->xsz + i + 1];
				NUMERIC_TYPE forward_slope = U.HV1x[j * Parptr->xsz + i + 1];

				NUMERIC_TYPE& local_slopelim =
					Arrptr->HV1x_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}
		}
	}

#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		for (int i = 0; i < Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z1x = Arrptr->DEM1x[j * Parptr->xsz + i];
			NUMERIC_TYPE ETA1x = Arrptr->ETA1x_slopelim[j * Parptr->xsz + i];

			NUMERIC_TYPE& H1x = U.H1x[j * Parptr->xsz + i];
			H1x = ETA1x - Z1x;
		}
	}

	const size_t size = Parptr->xsz * Parptr->ysz * sizeof(NUMERIC_TYPE);
	memcpy(U.HU1x, Arrptr->HU1x_slopelim, size);
	memcpy(U.HV1x, Arrptr->HV1x_slopelim, size);
}

void dg2::apply_slope_limiter_y
(
	Pars* Parptr,
	Solver* Solverptr,
	Arrays* Arrptr,
	FlowCoefficients const& U
)
{
#pragma omp parallel for
	for (int i = 0; i < Parptr->xsz; i++)
	{
		// north
		{
			const int j = 0;
			Arrptr->ETA1y_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HU1y_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HV1y_slopelim[j * Parptr->xsz + i] = C(0.0);
		}

		// south
		for (int i = 0; i < Parptr->xsz; i++)
		{
			const int j = Parptr->ysz - 1;
			Arrptr->ETA1y_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HU1y_slopelim[j * Parptr->xsz + i] = C(0.0);
			Arrptr->HV1y_slopelim[j * Parptr->xsz + i] = C(0.0);
		}
	}

#pragma omp parallel for
	for (int j = 1; j < Parptr->ysz - 1; j++)
	{
		for (int i = 0; i < Parptr->xsz; i++)
		{
			NUMERIC_TYPE minH = stencil_minH_y(Parptr, U, i, j);

			{
				NUMERIC_TYPE backward_const = eta(Parptr, Arrptr, U, i, j + 1);
				NUMERIC_TYPE backward_slope = eta1y(Parptr, Arrptr, U, i, j + 1);
				NUMERIC_TYPE local_const = eta(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE local_slope = eta1y(Parptr, Arrptr, U, i, j);
				NUMERIC_TYPE forward_const = eta(Parptr, Arrptr, U, i, j - 1);
				NUMERIC_TYPE forward_slope = eta1y(Parptr, Arrptr, U, i, j - 1);

				NUMERIC_TYPE& local_slopelim =
					Arrptr->ETA1y_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}

			{
				NUMERIC_TYPE backward_const = U.HU[(j + 1) * Parptr->xsz + i];
				NUMERIC_TYPE backward_slope = U.HU1y[(j + 1) * Parptr->xsz + i];
				NUMERIC_TYPE local_const = U.HU[j * Parptr->xsz + i];
				NUMERIC_TYPE local_slope = U.HU1y[j * Parptr->xsz + i];
				NUMERIC_TYPE forward_const = U.HU[(j - 1) * Parptr->xsz + i];
				NUMERIC_TYPE forward_slope = U.HU1y[(j - 1) * Parptr->xsz + i];

				NUMERIC_TYPE& local_slopelim =
					Arrptr->HU1y_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}

			{
				NUMERIC_TYPE backward_const = U.HV[(j + 1) * Parptr->xsz + i];
				NUMERIC_TYPE backward_slope = U.HV1y[(j + 1) * Parptr->xsz + i];
				NUMERIC_TYPE local_const = U.HV[j * Parptr->xsz + i];
				NUMERIC_TYPE local_slope = U.HV1y[j * Parptr->xsz + i];
				NUMERIC_TYPE forward_const = U.HV[(j - 1) * Parptr->xsz + i];
				NUMERIC_TYPE forward_slope = U.HV1y[(j - 1) * Parptr->xsz + i];

				NUMERIC_TYPE& local_slopelim =
					Arrptr->HV1y_slopelim[j * Parptr->xsz + i];
				local_slopelim = limit_slope(Solverptr, minH, Parptr->dx,
					backward_const, backward_slope,
					local_const, local_slope,
					forward_const, forward_slope);
			}
		}
	}

#pragma omp parallel for
	for (int j = 0; j < Parptr->ysz; j++)
	{
		for (int i = 0; i < Parptr->xsz; i++)
		{
			NUMERIC_TYPE Z1y = Arrptr->DEM1y[j * Parptr->xsz + i];
			NUMERIC_TYPE ETA1y = Arrptr->ETA1y_slopelim[j * Parptr->xsz + i];

			NUMERIC_TYPE& H1y = U.H1y[j * Parptr->xsz + i];
			H1y = ETA1y - Z1y;
		}
	}

	const size_t size = Parptr->xsz * Parptr->ysz * sizeof(NUMERIC_TYPE);
	memcpy(U.HU1y, Arrptr->HU1y_slopelim, size);
	memcpy(U.HV1y, Arrptr->HV1y_slopelim, size);
}

NUMERIC_TYPE dg2::limit_slope
(
	Solver* Solverptr,
	NUMERIC_TYPE stencil_minH,
	NUMERIC_TYPE mesh_delta,
	NUMERIC_TYPE backward_const,
	NUMERIC_TYPE backward_slope,
	NUMERIC_TYPE local_const,
	NUMERIC_TYPE local_slope,
	NUMERIC_TYPE forward_const,
	NUMERIC_TYPE forward_slope
)
{
	if (stencil_minH <= C(0.30) * Solverptr->maxH) return local_slope;

	NUMERIC_TYPE backward_neg = limit_neg(backward_const, backward_slope);
	NUMERIC_TYPE local_pos = limit_pos(local_const, local_slope);
	NUMERIC_TYPE local_neg = limit_neg(local_const, local_slope);
	NUMERIC_TYPE forward_pos = limit_pos(forward_const, forward_slope);

	NUMERIC_TYPE backward_jump = FABS(local_pos - backward_neg);
	NUMERIC_TYPE forward_jump = FABS(local_neg - forward_pos);

	NUMERIC_TYPE norm = FMAX(
		FABS(local_const - backward_const), FABS(forward_const - local_const));

	NUMERIC_TYPE DS_backward = C(0.0);
	NUMERIC_TYPE DS_forward = C(0.0);

	if (FABS(norm) > C(1e-12))
	{
		DS_backward = backward_jump / (C(0.5) * mesh_delta * norm);
		DS_forward = forward_jump / (C(0.5) * mesh_delta * norm);
	}

	if (FMAX(DS_backward, DS_forward) >= Solverptr->krivodonova_threshold)
	{
		return minmod
		(
			local_slope,
			(local_const - backward_const) / SQRT(C(3.0)),
			(forward_const - local_const) / SQRT(C(3.0))
		);
	}
	else
	{
		return local_slope;
	}
}

NUMERIC_TYPE dg2::minmod(NUMERIC_TYPE a, NUMERIC_TYPE b, NUMERIC_TYPE c)
{
	if (a * b > C(0.0) && a * c > C(0.0))
	{
		NUMERIC_TYPE sign_a = (a > C(0.0)) - (a < C(0.0));
		return sign_a * (std::min)({ FABS(a), FABS(b), FABS(c) });
	}
	else
	{
		return C(0.0);
	}
}

NUMERIC_TYPE dg2::stencil_minH_x
(
	Pars* Parptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE backward_const = U.H[j * Parptr->xsz + i - 1];
	NUMERIC_TYPE backward_slope = U.H1x[j * Parptr->xsz + i - 1];
	NUMERIC_TYPE local_const = U.H[j * Parptr->xsz + i];
	NUMERIC_TYPE local_slope = U.H1x[j * Parptr->xsz + i];
	NUMERIC_TYPE forward_const = U.H[j * Parptr->xsz + i + 1];
	NUMERIC_TYPE forward_slope = U.H1x[j * Parptr->xsz + i + 1];
	return (std::min)({
				gauss_lower(backward_const, backward_slope),
				gauss_upper(backward_const, backward_slope),
				gauss_lower(local_const, local_slope),
				gauss_upper(local_const, local_slope),
				gauss_lower(forward_const, forward_slope),
				gauss_upper(forward_const, forward_slope)
		});
}

NUMERIC_TYPE dg2::stencil_minH_y
(
	Pars* Parptr,
	FlowCoefficients const& U,
	const int i,
	const int j
)
{
	NUMERIC_TYPE backward_H = U.H[(j + 1) * Parptr->xsz + i];
	NUMERIC_TYPE backward_H1x = U.H1x[(j + 1) * Parptr->xsz + i];
	NUMERIC_TYPE local_H = U.H[j * Parptr->xsz + i];
	NUMERIC_TYPE local_H1x = U.H1x[j * Parptr->xsz + i];
	NUMERIC_TYPE forward_H = U.H[(j - 1) * Parptr->xsz + i];
	NUMERIC_TYPE forward_H1x = U.H1x[(j - 1) * Parptr->xsz + i];

	return (std::min)({
			gauss_lower(backward_H, backward_H1x),
			gauss_upper(backward_H, backward_H1x),
			gauss_lower(local_H, local_H1x),
			gauss_upper(local_H, local_H1x),
			gauss_lower(forward_H, forward_H1x),
			gauss_upper(forward_H, forward_H1x)
		});
}
