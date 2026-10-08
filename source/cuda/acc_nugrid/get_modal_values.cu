#include "get_modal_values.cuh"

__host__ void lis::cuda::acc_nugrid::get_modal_values
(
	AssembledSolution& d_assem_sol,
	const int&         mesh_dim,
	const int&         interface_dim,
	const Fnames&      filenames,
	const States&      states,
	const NUMERIC_TYPE& no_data
)
{

//		char qx_raster_filename[64];
//		char qy_raster_filename[64];
		char dem_filename[64];
		char dem1x_filename[64];
		char dem1y_filename[64];
		
//		strcpy(qx_raster_filename, start_file);
//		strcpy(qy_raster_filename, start_file);
		strcpy(dem_filename, filenames.demfilename);
		strcpy(dem1x_filename, filenames.demfilename);
		strcpy(dem1y_filename, filenames.demfilename);
		
		//strcat(h_raster_filename,  ".start");
//		strcat(qx_raster_filename, ".Qx");
//		strcat(qy_raster_filename, ".Qy");
		//strcat(dem_filename,       ".dem");
		strcat(dem1x_filename, "1x");
		strcat(dem1y_filename, "1y");

		NUMERIC_TYPE* h_raster  = new NUMERIC_TYPE[d_assem_sol.length]();
		NUMERIC_TYPE* v_raster = new NUMERIC_TYPE[d_assem_sol.length]();
//		NUMERIC_TYPE* q1_raster = new NUMERIC_TYPE[d_assem_sol.length]();
//		NUMERIC_TYPE* q2_raster = new NUMERIC_TYPE[d_assem_sol.length]();
//		NUMERIC_TYPE* q3_raster = new NUMERIC_TYPE[d_assem_sol.length]();
//		NUMERIC_TYPE* q4_raster = new NUMERIC_TYPE[d_assem_sol.length]();
		NUMERIC_TYPE* dem       = new NUMERIC_TYPE[d_assem_sol.length]();
		NUMERIC_TYPE* dem1x     = new NUMERIC_TYPE[d_assem_sol.length]();
		NUMERIC_TYPE* dem1y     = new NUMERIC_TYPE[d_assem_sol.length]();
//		NUMERIC_TYPE* demxy     = new NUMERIC_TYPE[d_assem_sol.length]();
		
		for (int i = 0; i < d_assem_sol.length; i++) dem[i] = no_data; // C(100.0); // high wall, diff set up for open BCs
//		for (int i = 0; i < d_assem_sol.length; i++) dem1x[i] = C(0.0); // high wall, diff set up for open BCs // No need cause it's zero by default
//		for (int i = 0; i < d_assem_sol.length; i++) dem1y[i] = C(0.0); // high wall, diff set up for open BCs // No need cause it's zero by default
	
//		read_raster_file(qx_raster_filename, q1_raster, mesh_dim);
//		read_raster_file(qy_raster_filename, q2_raster, mesh_dim);
//		read_raster_file(qx_raster_filename, q3_raster, mesh_dim);
//		read_raster_file(qy_raster_filename, q4_raster, mesh_dim);
		read_raster_file(dem_filename,       dem,       mesh_dim, no_data);
		read_raster_file(dem1x_filename,     dem1x,     mesh_dim, no_data);
		read_raster_file(dem1y_filename,     dem1y,     mesh_dim, no_data);
		
		if (states.startfile == ON)
		{
			char h_raster_filename[64];

			strcpy(h_raster_filename, filenames.startfilename);

			read_raster_file(h_raster_filename, h_raster, mesh_dim, C(0.0));

			//if (states.startelev == ON)
			//{
			//	StartFile::subtract_dem(H, DEM, geometry, pitch, offset);
			//}
		}


		//int nrows = 0;
		//int ncols = 0;

		//NUMERIC_TYPE dummy = C(0.0);

		//FILE* fp = fopen(dem_filename, "r");

		//if (NULL == fp)
		//{
		//	fprintf(stdout, "No raster file found: %s, using default values.\n", dem_filename);
		//	return;
		//}

		//fscanf(fp, "%s %d", buf, &ncols);
		//fscanf(fp, "%s %d", buf, &nrows);
		//fscanf(fp, "%s %" NUM_FRMT, buf, &dummy);
		//fscanf(fp, "%s %" NUM_FRMT, buf, &dummy);
		//fscanf(fp, "%s %" NUM_FRMT, buf, &dummy);
		//fscanf(fp, "%s %" NUM_FRMT, buf, &dummy);

		//fclose(fp);

		size_t bytes = d_assem_sol.length * sizeof(NUMERIC_TYPE);

		copy_cuda(d_assem_sol.h,   h_raster, bytes);
		copy_cuda(d_assem_sol.v,   v_raster, bytes);
//		copy_cuda(d_assem_sol.q1,  q1_raster, bytes);
//		copy_cuda(d_assem_sol.q2,  q2_raster, bytes);
//		copy_cuda(d_assem_sol.q3,  q3_raster, bytes);
//		copy_cuda(d_assem_sol.q4,  q4_raster, bytes);
		copy_cuda(d_assem_sol.z0,  dem,   bytes);
		copy_cuda(d_assem_sol.z1x, dem1x, bytes);
		copy_cuda(d_assem_sol.z1y, dem1y, bytes);
//		copy_cuda(d_assem_sol.zxy, demxy,     bytes);
		
		if (strlen(filenames.nfilename) > 0) {
			char n_filename[64];
			strcpy(n_filename, filenames.nfilename);
			NUMERIC_TYPE* n_raster = new NUMERIC_TYPE[d_assem_sol.length]();
			for (int i = 0; i < d_assem_sol.length; i++) n_raster[i] = C(0.0); // high wall, diff set up for open BCs
			read_raster_file(n_filename, n_raster, mesh_dim, C(0.0));
			copy_cuda(d_assem_sol.n0, n_raster, bytes);
			delete[] n_raster;
		}

		delete[] h_raster;
		delete[] v_raster;
//		delete[] q1_raster;
//		delete[] q2_raster;
//		delete[] q3_raster;
//		delete[] q4_raster;
		delete[] dem;
		delete[] dem1x;
		delete[] dem1y;
//		delete[] demxy;

//	}
}