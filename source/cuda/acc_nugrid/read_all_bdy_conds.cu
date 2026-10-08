#include "read_all_bdy_conds.cuh"

void lis::cuda::acc_nugrid::read_all_bdy_conds
(
	const char* bcifilename,
	Boundaries& boundaries,
	const Pars& pars
)
{
	read_bdy_conds(bcifilename, NORTH, boundaries.north, pars);
	read_bdy_conds(bcifilename, EAST, boundaries.east, pars);
	read_bdy_conds(bcifilename, SOUTH, boundaries.south, pars);
	read_bdy_conds(bcifilename, WEST, boundaries.west, pars);
}