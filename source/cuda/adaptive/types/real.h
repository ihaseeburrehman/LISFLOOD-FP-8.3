#pragma once

#if _NUMERIC_MODE == 1

typedef double real;
#define C(x) x
#define NUM_FRMT "lf"
#define NUM_FIG  ".15"

#else

typedef float real;
#define C(x) x##f
#define NUM_FRMT "f"
#define NUM_FIG  ".8"

#endif