#pragma once

#include "../lisflood.h"

void read_ascfile
(
	const char* filename,
	Pars *Parptr,
	NUMERIC_TYPE *field,
	const char* message,
	const int verbose
);

