#ifndef LIS_WINDOWS_EXT_H_
#define LIS_WINDOWS_EXT_H_

#include <windows.h> 

struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz);

#endif