#pragma once

#include <time.h>

#if defined(_MSC_VER)
#include "windows_ext.h"
#else
#include <sys/time.h>
#endif

#define GRID_ALIGN_WIDTH 8
#define	THREAD_SG_CELL_BLOCK 64