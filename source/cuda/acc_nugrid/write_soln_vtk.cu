#include "write_soln_vtk.cuh"

__host__ void lis::cuda::acc_nugrid::write_soln_vtk
(
	const char*                 respath,
	const AssembledSolution&    d_assem_sol,
	const Pars& pars,
	const int&                  lev,
	const NUMERIC_TYPE&                 depth_thresh,
	const int&                  call_gzip,
	const int& precision
)
{
	char fullpath[800] = {'\0'};
	char tmp_sys_com[255];

	
	sprintf(fullpath, "%s-%.4d%s", respath, pars.SaveNo , ".vtk");

	FILE* fp = fopen(fullpath, "w");

	if (NULL == fp)
	{
		fprintf(stderr, "Error writing VTK file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	fprintf
	(
		fp,
		"# vtk DataFile Version 3.0\n"
		"Multiresolution flow and topo data\n"
		"ASCII\n"
		"\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS %d float\n",
		d_assem_sol.length * 4
	);

	// xmax, ymax from params may have non-zero origin
	// however, Morton codes assume zero origin
	// hence, modified xmax, ymax for bounds checking
//	NUMERIC_TYPE xmax_0_orig = params.xmin + params.xsz * dx_finest;
//	NUMERIC_TYPE ymax_0_orig = params.ymin + params.ysz * dx_finest;

	NUMERIC_TYPE xmax_0_orig = pars.xsz * pars.dx;
	NUMERIC_TYPE ymax_0_orig = pars.ysz * pars.dx;

	NUMERIC_TYPE*     z0        = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE*     h         = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE*     v         = new NUMERIC_TYPE[d_assem_sol.length];
	index_1D* act_idcs = new index_1D[d_assem_sol.length];
	int*      levels   = new int[d_assem_sol.length];

	size_t bytes_flow     = d_assem_sol.length * sizeof(NUMERIC_TYPE);
	size_t bytes_act_idcs = d_assem_sol.length * sizeof(index_1D);
	size_t bytes_levels   = d_assem_sol.length * sizeof(int);
	
	copy_cuda(z0,        d_assem_sol.z0, bytes_flow);
	copy_cuda(h,         d_assem_sol.h, bytes_flow);
	copy_cuda(v,        d_assem_sol.v, bytes_flow);
	copy_cuda(act_idcs, d_assem_sol.act_idcs, bytes_act_idcs);
	copy_cuda(levels,   d_assem_sol.levels, bytes_levels);

	// number of cells excluding those in extended domain
	int num_bound = 0;

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		index_1D act_idx = act_idcs[i];
		int      level   = levels[i];

		MortonCode code = act_idx - get_lvl_idx(level);

		NUMERIC_TYPE local_cell_size = pars.dx * ( 1 << (lev - level) );

		int x = compact(code);
		int y = compact(code >> 1);

//		NUMERIC_TYPE x_centre = params.xmin + (x * local_cell_size + local_cell_size / C(2.0));
//		NUMERIC_TYPE y_centre = params.ymax - (y * local_cell_size + local_cell_size / C(2.0));

		NUMERIC_TYPE x_centre = x * local_cell_size + local_cell_size / C(2.0);
		NUMERIC_TYPE y_centre = y * local_cell_size + local_cell_size / C(2.0);

		bool bound = (x_centre < xmax_0_orig && y_centre < ymax_0_orig);

		//if (!bound) continue;

		num_bound++;

		Points points =
		{
			pars.blx + x * local_cell_size,       // lower left  x
			pars.bly + y * local_cell_size,       // lower left  y
			pars.blx + x * local_cell_size,       // upper left  x
			pars.bly + (y + 1) * local_cell_size, // upper left  y
			pars.blx + (x + 1) * local_cell_size, // lower right x
			pars.bly + y * local_cell_size,       // lower right y
			pars.blx + (x + 1) * local_cell_size, // upper right x
			pars.bly + (y + 1) * local_cell_size  // upper right y

//			params.xmin + x * local_cell_size,       // lower left  x
//			params.ymax - y * local_cell_size,       // lower left  y
//			params.xmin + x * local_cell_size,       // upper left  x
//			params.ymax - (y + 1) * local_cell_size, // upper left  y
//			params.xmin + (x + 1) * local_cell_size, // lower right x
//			params.ymax - y * local_cell_size,       // lower right y
//			params.xmin + (x + 1) * local_cell_size, // upper right x
//			params.ymax - (y + 1) * local_cell_size  // upper right y
		};

		fprintf
		(
			fp,
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n",
			precision, points.ll_x, precision, points.ll_y, precision, C(1.0),
			precision, points.ul_x, precision, points.ul_y, precision, C(1.0),
			precision, points.lr_x, precision, points.lr_y, precision, C(1.0),
			precision, points.ur_x, precision, points.ur_y, precision, C(1.0)
		);
	}

	fprintf(fp, "\nCELLS %d %d\n", d_assem_sol.length, d_assem_sol.length * 5);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		// point counter to make sure correct vertices per cell
		int pt_ctr = i * 4;

		fprintf(fp, "4 %d %d %d %d\n", 0 + pt_ctr, 1 + pt_ctr, 2 + pt_ctr, 3 + pt_ctr);
	}

	fprintf(fp, "\nCELL_TYPES %d\n", d_assem_sol.length);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		fprintf(fp, "%d\n", 8);
	}

	fprintf(fp, "\nCELL_DATA %d\n", d_assem_sol.length);

	fprintf
	(
		fp,
		"\nSCALARS z0 float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
//		fprintf( fp, "%" NUM_FMT "\n", z0[i] );
		fprintf(fp, "%.*" NUM_FMT "\n", precision, z0[i]);
	}

	fprintf
	(
		fp,
		"\nSCALARS h float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
//		fprintf(fp, "%" NUM_FMT "\n", h[i]);

		fprintf(fp, "%.*" NUM_FMT "\n", precision, h[i]);
	}

	fprintf
	(
		fp,
		"\nSCALARS v float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		//		fprintf(fp, "%" NUM_FMT "\n", h[i]);

		fprintf(fp, "%.*" NUM_FMT "\n", precision, v[i]);
	}

	//for (int i = 0; i < d_assem_sol.length; i++)
	//{
	//	NUMERIC_TYPE vx = (h[i] < depth_thresh) ? C(0.0) : qx[i] / h[i];
	//	
	//	fprintf( fp, "%" NUM_FRMT "\n", vx );
	//}

	//fprintf
	//(
	//	fp,
	//	"\nSCALARS vy float 1\n"
	//	"LOOKUP_TABLE default\n"
	//);

	//for (int i = 0; i < d_assem_sol.length; i++)
	//{
	//	NUMERIC_TYPE vy = (h[i] < depth_thresh) ? C(0.0) : qy[i] / h[i];
	//	
	//	fprintf( fp, "%" NUM_FRMT "\n", vy );
	//}

	fclose(fp);

	// check if we need to zip the file up
	if (call_gzip == ON)
	{
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", fullpath);
		system(tmp_sys_com);
	}


	delete[] z0;
	delete[] h;
	delete[] v;
	delete[] act_idcs;
	delete[] levels;
}

__host__ void lis::cuda::acc_nugrid::write_soln_vtk_with_n
(
	const char* respath,
	const AssembledSolution& d_assem_sol,
	const Pars& pars,
	const int& lev,
	const NUMERIC_TYPE& depth_thresh,
	const int& call_gzip,
	const int& precision
)
{
	char fullpath[255] = { '\0' };
	char tmp_sys_com[255];

	sprintf(fullpath, "%s-%.4d%s", respath, pars.SaveNo, ".vtk");

	FILE* fp = fopen(fullpath, "w");

	if (NULL == fp)
	{
		fprintf(stderr, "Error writing VTK file, file: %s, line: %d.\n", __FILE__, __LINE__);
		exit(-1);
	}

	fprintf
	(
		fp,
		"# vtk DataFile Version 3.0\n"
		"Multiresolution flow and topo data\n"
		"ASCII\n"
		"\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS %d float\n",
		d_assem_sol.length * 4
	);

	// xmax, ymax from params may have non-zero origin
	// however, Morton codes assume zero origin
	// hence, modified xmax, ymax for bounds checking
//	NUMERIC_TYPE xmax_0_orig = params.xmin + params.xsz * dx_finest;
//	NUMERIC_TYPE ymax_0_orig = params.ymin + params.ysz * dx_finest;

	NUMERIC_TYPE xmax_0_orig = pars.xsz * pars.dx;
	NUMERIC_TYPE ymax_0_orig = pars.ysz * pars.dx;

	NUMERIC_TYPE* z0 = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE* h = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE* v = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE* n0 = new NUMERIC_TYPE[d_assem_sol.length];
	index_1D* act_idcs = new index_1D[d_assem_sol.length];
	int* levels = new int[d_assem_sol.length];

	size_t bytes_flow = d_assem_sol.length * sizeof(NUMERIC_TYPE);
	size_t bytes_act_idcs = d_assem_sol.length * sizeof(index_1D);
	size_t bytes_levels = d_assem_sol.length * sizeof(int);

	copy_cuda(z0, d_assem_sol.z0, bytes_flow);
	copy_cuda(h, d_assem_sol.h, bytes_flow);
	copy_cuda(v, d_assem_sol.v, bytes_flow);
	copy_cuda(n0, d_assem_sol.n0, bytes_flow);
	copy_cuda(act_idcs, d_assem_sol.act_idcs, bytes_act_idcs);
	copy_cuda(levels, d_assem_sol.levels, bytes_levels);

	// number of cells excluding those in extended domain
	int num_bound = 0;

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		index_1D act_idx = act_idcs[i];
		int      level = levels[i];

		MortonCode code = act_idx - get_lvl_idx(level);

		NUMERIC_TYPE local_cell_size = pars.dx * (1 << (lev - level));

		int x = compact(code);
		int y = compact(code >> 1);

		//		NUMERIC_TYPE x_centre = params.xmin + (x * local_cell_size + local_cell_size / C(2.0));
		//		NUMERIC_TYPE y_centre = params.ymax - (y * local_cell_size + local_cell_size / C(2.0));

		NUMERIC_TYPE x_centre = x * local_cell_size + local_cell_size / C(2.0);
		NUMERIC_TYPE y_centre = y * local_cell_size + local_cell_size / C(2.0);

		bool bound = (x_centre < xmax_0_orig&& y_centre < ymax_0_orig);

		//if (!bound) continue;

		num_bound++;

		Points points =
		{
			pars.blx + x * local_cell_size,       // lower left  x
			pars.bly + y * local_cell_size,       // lower left  y
			pars.blx + x * local_cell_size,       // upper left  x
			pars.bly + (y + 1) * local_cell_size, // upper left  y
			pars.blx + (x + 1) * local_cell_size, // lower right x
			pars.bly + y * local_cell_size,       // lower right y
			pars.blx + (x + 1) * local_cell_size, // upper right x
			pars.bly + (y + 1) * local_cell_size  // upper right y

//			params.xmin + x * local_cell_size,       // lower left  x
//			params.ymax - y * local_cell_size,       // lower left  y
//			params.xmin + x * local_cell_size,       // upper left  x
//			params.ymax - (y + 1) * local_cell_size, // upper left  y
//			params.xmin + (x + 1) * local_cell_size, // lower right x
//			params.ymax - y * local_cell_size,       // lower right y
//			params.xmin + (x + 1) * local_cell_size, // upper right x
//			params.ymax - (y + 1) * local_cell_size  // upper right y
		};

		fprintf
		(
			fp,
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n"
			"%.*" NUM_FMT " %.*" NUM_FMT " %.*" NUM_FMT "\n",
			precision, points.ll_x, precision, points.ll_y, precision, C(1.0),
			precision, points.ul_x, precision, points.ul_y, precision, C(1.0),
			precision, points.lr_x, precision, points.lr_y, precision, C(1.0),
			precision, points.ur_x, precision, points.ur_y, precision, C(1.0)
		);
	}

	fprintf(fp, "\nCELLS %d %d\n", d_assem_sol.length, d_assem_sol.length * 5);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		// point counter to make sure correct vertices per cell
		int pt_ctr = i * 4;

		fprintf(fp, "4 %d %d %d %d\n", 0 + pt_ctr, 1 + pt_ctr, 2 + pt_ctr, 3 + pt_ctr);
	}

	fprintf(fp, "\nCELL_TYPES %d\n", d_assem_sol.length);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		fprintf(fp, "%d\n", 8);
	}

	fprintf(fp, "\nCELL_DATA %d\n", d_assem_sol.length);

	fprintf
	(
		fp,
		"\nSCALARS z0 float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
//		fprintf(fp, "%" NUM_FMT "\n", z0[i]);
		fprintf(fp, "%.*" NUM_FMT "\n", precision, z0[i]);
	}

	fprintf
	(
		fp,
		"\nSCALARS h float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
//		fprintf(fp, "%" NUM_FMT "\n", h[i]);
		fprintf(fp, "%.*" NUM_FMT "\n", precision, h[i]);
	}

	fprintf
	(
		fp,
		"\nSCALARS v float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
		//		fprintf(fp, "%" NUM_FMT "\n", h[i]);
		fprintf(fp, "%.*" NUM_FMT "\n", precision, v[i]);
	}

	fprintf
	(
		fp,
		"\nSCALARS n float 1\n"
		"LOOKUP_TABLE default\n"
	);

	for (int i = 0; i < d_assem_sol.length; i++)
	{
//		fprintf(fp, "%" NUM_FMT "\n", n0[i]);
		fprintf(fp, "%.*" NUM_FMT "\n", precision, n0[i]);
	}

	//fprintf
	//(
	//	fp,
	//	"\nSCALARS vx float 1\n"
	//	"LOOKUP_TABLE default\n"
	//);

	//for (int i = 0; i < d_assem_sol.length; i++)
	//{
	//	NUMERIC_TYPE vx = (h[i] < depth_thresh) ? C(0.0) : qx[i] / h[i];
	//	
	//	fprintf( fp, "%" NUM_FRMT "\n", vx );
	//}

	//fprintf
	//(
	//	fp,
	//	"\nSCALARS vy float 1\n"
	//	"LOOKUP_TABLE default\n"
	//);

	//for (int i = 0; i < d_assem_sol.length; i++)
	//{
	//	NUMERIC_TYPE vy = (h[i] < depth_thresh) ? C(0.0) : qy[i] / h[i];
	//	
	//	fprintf( fp, "%" NUM_FRMT "\n", vy );
	//}

	fclose(fp);

	// check if we need to zip the file up
	if (call_gzip == ON)
	{
		sprintf(tmp_sys_com, "%s%s", "gzip -9 -f ", fullpath);
		system(tmp_sys_com);
	}

	delete[] z0;
	delete[] h;
	delete[] v;
	delete[] n0;
	delete[] act_idcs;
	delete[] levels;
}