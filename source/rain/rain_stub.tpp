#include "rain.h"
#include "../utility.h"

template<class Allocator>
DynamicRain<Allocator>::DynamicRain
(
    const char* filename,
    int verbose,
    const Allocator& allocator
)
{}

template<class Allocator>
bool DynamicRain<Allocator>::has_same_origin
(
    Pars* Parptr
)
{
    return false;
}

template<class Allocator>
bool DynamicRain<Allocator>::is_tile_size_multiple_of_grid
(
    Pars* Parptr
)
{
    return false;
}

template<class Allocator>
void DynamicRain<Allocator>::update_time
(
    NUMERIC_TYPE t
)
{}

template<class Allocator>
bool DynamicRain<Allocator>::update_time_SGC
(
    NUMERIC_TYPE t
)
{
    return false;
}

template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell
(
    Pars* Parptr,
    int i,
    int j
)
{
    throw std::logic_error("Not implemented");
}

template<class Allocator>
NUMERIC_TYPE DynamicRain<Allocator>::rate_at_cell_SGC
(
    int i,
    int j,
    const NUMERIC_TYPE dx_col, 
    const NUMERIC_TYPE dy_col,
    NUMERIC_TYPE tly
)
{
    throw std::logic_error("Not implemented");
}

template<class Allocator>
void DynamicRain<Allocator>::update_H
(
    Pars* Parptr,
    Solver *Solverptr,
    Arrays *Arrptr
)
{}

template<class Allocator>
void DynamicRain<Allocator>::update_rain_grid_SGM
(
    const NUMERIC_TYPE curr_time,
    NUMERIC_TYPE* rain_grid,
    const NUMERIC_TYPE* dem_grid,
    WetDryRowBound* wet_dry_bounds,
    const NUMERIC_TYPE* dx_col,
    const NUMERIC_TYPE* dy_col,
    const NUMERIC_TYPE tyl,
    const int grid_rows,
    const int grid_cols_padded,
    const NUMERIC_TYPE* cell_area_col
)
{}

template<class Allocator>
DynamicRain<Allocator>::~DynamicRain()
{}
