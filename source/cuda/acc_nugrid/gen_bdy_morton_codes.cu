#include "gen_bdy_morton_codes.cuh"

void lis::cuda::acc_nugrid::gen_bdy_morton_codes
(
	const Boundary& boundary,
	const Pars& pars,
	const int& direction
)
{
	int current = 0;
	
	switch (direction)
	{
	case SOUTH:
	{
		for (int i = 0; i < boundary.num_cells(); i++)
		{
			current = boundary.start + i;

//			boundary.codes[i] = generate_morton_code(current, params.ysz - 1);
			boundary.codes[i] = generate_morton_code(current, 0);
		}

		break;
	}
	case NORTH:
	{
		for (int i = 0; i < boundary.num_cells(); i++)
		{
			current = boundary.start + i;

//			boundary.codes[i] = generate_morton_code(current, 0);
			boundary.codes[i] = generate_morton_code(current, pars.ysz - 1);
		}

		break;
	}
	case EAST:
	{
		for (int i = 0; i < boundary.num_cells(); i++)
		{
			current = boundary.start + i;

			boundary.codes[i] = generate_morton_code(pars.xsz - 1, current);
		}

		break;
	}
	case WEST:
	{
		for (int i = 0; i < boundary.num_cells(); i++)
		{
			current = boundary.start + i;

			boundary.codes[i] = generate_morton_code(0, current);
		}

		break;
	}
	default:
		break;
	}
}