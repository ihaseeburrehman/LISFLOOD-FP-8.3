#include "read_bdy_conds.cuh"

void lis::cuda::acc_nugrid::read_bdy_conds
(
	const char*                 bcifilename,
	const int                   direction,
	Boundary&                   boundary,
	const Pars&                 pars
)
{
	NUMERIC_TYPE upper = C(0.0);
	NUMERIC_TYPE lower = C(0.0);

	char bcidir  = '0';
	char filedir = '0';

	char str[255];
	char bdytype_buf[8];

	NUMERIC_TYPE origin = C(0.0);

//	int BCi = 0;

	switch (direction)
	{
	case NORTH:
		origin = pars.blx;
		bcidir = 'N';
		break;
	case EAST:
		origin = pars.bly;
		bcidir = 'E';
		break;
	case SOUTH:
		origin = pars.blx;
		bcidir = 'S';
		break;
	case WEST:
		origin = pars.bly;
		bcidir = 'W';
		break;
	default:
		break;
	}
	
	FILE* fp = fopen(bcifilename, "r");

	if (NULL == fp)
	{
		fprintf(stderr, "Error opening boundary condition file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	while (filedir != bcidir)
	{
		if ( NULL == fgets(str, sizeof(str), fp) )
		{
			fprintf(stdout, "No specifications found for boundary %c, proceeding with closed boundaries.\n", bcidir);
			fclose(fp);

			return;
		}

		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s", &filedir, &lower, &upper, bdytype_buf);
	}

	fclose(fp);

	if (bcidir == 'N') {
//		int BCi1 = -1;

		if (lower < pars.blx) lower = pars.blx;
		if (lower > pars.blx + pars.xsz * pars.dx) lower = pars.blx + pars.xsz * pars.dx;
		if (upper < pars.blx) upper = pars.blx;
		if (upper > pars.blx + pars.xsz * pars.dx) upper = pars.blx + pars.xsz * pars.dx;

//		BCi1 = (int)((lower - pars.blx) / pars.dy);

//		BCi = BCi1;
	} 

	if (bcidir == 'W') {

//		int BCi1 = -1;

		if (lower < pars.bly) lower = pars.bly;
		if (lower > pars.tly) lower = pars.tly;
		if (upper < pars.bly) upper = pars.bly;
		if (upper > pars.tly) upper = pars.tly;

//		BCi1 = (int)(2 * pars.xsz + 2 * pars.ysz - (pars.tly - lower) / pars.dy);

//		BCi = BCi1;
	}
	
	if (bcidir == 'S')
	{
//		int BCi1 = -1;

		if (lower < pars.blx) lower = pars.blx;
		if (lower > pars.blx + pars.xsz * pars.dx) lower = pars.blx + pars.xsz * pars.dx;
		if (upper < pars.blx) upper = pars.blx;
		if (upper > pars.blx + pars.xsz * pars.dx) upper = pars.blx + pars.xsz * pars.dx;

//		BCi1 = (int)(2 * pars.xsz + pars.ysz - (lower - pars.blx) / pars.dy);

//		BCi = BCi1;
	}

	if (bcidir == 'E')
	{
//		int BCi1 = -1;

		if (lower < pars.bly) lower = pars.bly;
		if (lower > pars.tly) lower = pars.tly;
		if (upper < pars.bly) upper = pars.bly;
		if (upper > pars.tly) upper = pars.tly;

//		BCi1 = (int)(pars.xsz + (pars.tly - lower) / pars.dy);

//		BCi = BCi1;
	}


	int  bdytype        = CLOSED;
	NUMERIC_TYPE inlet          = C(0.0);
	
	char timeseries[32] = {'\0'};

	if ( !strncmp(bdytype_buf, "CLOSED", 6) )
	{
		bdytype = CLOSED;
	}
	else if ( !strncmp(bdytype_buf, "FREE", 4) )
	{
		bdytype = FREE;
	}
	else if ( !strncmp(bdytype_buf, "HFIX", 4) )
	{
		bdytype = HFIX;
		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &filedir, &lower, &upper, bdytype_buf, &inlet);
	}
	else if ( !strncmp(bdytype_buf, "HVAR", 4) )
	{
		bdytype = HVAR;
		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &filedir, &lower, &upper, bdytype_buf, timeseries);
	}
	else if ( !strncmp(bdytype_buf, "QFIX", 4) )
	{
		bdytype = QFIX;
		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %" NUM_FMT, &filedir, &lower, &upper, bdytype_buf, &inlet);
	}
	else if ( !strncmp(bdytype_buf, "QVAR", 4) )
	{
		bdytype = QVAR;
		sscanf(str, "%c %" NUM_FMT " %" NUM_FMT " %s %s", &filedir, &lower, &upper, bdytype_buf, timeseries);
	}

	boundary.start = (lower - origin) / pars.dx;
	boundary.end = (upper - origin) / pars.dx - 1;
	gen_bdy_morton_codes(boundary, pars, direction);
	boundary.bdytype = bdytype;
	boundary.inlet = inlet;
//	boundary.BCi = BCi;
	sprintf(boundary.timeseries, "%s", timeseries);
}