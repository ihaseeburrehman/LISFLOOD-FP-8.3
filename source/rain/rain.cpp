#include "rain.h"

#if _NETCDF == 1
    #include "rain.tpp"
#else
    #include "rain_stub.tpp"
#endif

template class DynamicRain<>;
