#pragma once
#include "../lisflood.h"

NUMERIC_TYPE linear_interpolate
(
	TimeSeries *timeSeries,
	NUMERIC_TYPE t
);

int boundary_index_w
(
	Pars *Parptr,
	const int i,
	const int j
);

int boundary_index_e
(
	Pars *Parptr,
	const int i,
	const int j
);

int boundary_index_n
(
	Pars *Parptr,
	const int i,
	const int j
);

int boundary_index_s
(
	Pars *Parptr,
	const int i,
	const int j
);
