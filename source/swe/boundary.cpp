#include "boundary.h"

NUMERIC_TYPE linear_interpolate
(
	TimeSeries * timeSeries,
	NUMERIC_TYPE t
)
{
	const NUMERIC_TYPE * times = timeSeries->time;
	const NUMERIC_TYPE * values = timeSeries->value;

	if (t < times[0]) return values[0];

	for (int i=1; i < timeSeries->count; i++)
	{
		if (times[i - 1] <= t && times[i] > t)
		{
			NUMERIC_TYPE dt, a1, a2;
			dt = times[i] - times[i - 1];
			a1 = (t - times[i - 1]) / dt;
			a2 = C(1.0) - a1;

			return a1*values[i] + a2*values[i - 1];
		}
	}
	
	return values[timeSeries->count - 1];
}

int boundary_index_w
(
	Pars *Parptr,
	const int i,
	const int j
)
{
	return 2*Parptr->xsz + 2*Parptr->ysz - j - 1;
}

int boundary_index_e
(
	Pars *Parptr,
	const int i,
	const int j
)
{
	return Parptr->xsz + j;
}

int boundary_index_n
(
	Pars *Parptr,
	const int i,
	const int j
)
{
	return i;
}

int boundary_index_s
(
	Pars *Parptr,
	const int i,
	const int j
)
{
	return 2*Parptr->xsz + Parptr->ysz - i - 1;
}
