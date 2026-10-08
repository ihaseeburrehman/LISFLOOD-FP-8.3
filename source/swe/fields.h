#pragma once
#include "../lisflood.h"

void allocate_swe_fields
(
	Pars *Parptr,
	Arrays *Arrptr
);

void deallocate_swe_fields(Arrays *Arrptr);
