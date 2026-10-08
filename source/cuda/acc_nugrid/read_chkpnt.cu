#include "read_chkpnt.cuh"
#include "../../VersionHistory.h"

__host__ void lis::cuda::acc_nugrid::read_chkpnt
(
	AssembledSolution&    d_assem_sol,
	NonUniformNeighbours& d_non_uniform_nghbrs,
	Fnames& filenames,
	States& states,
	Pars& pars,
	Solver& solver,
	NUMERIC_TYPE* d_maxH,
	NUMERIC_TYPE* d_totalHtm,
	NUMERIC_TYPE* d_maxHtm,
	NUMERIC_TYPE* d_initHtm,
	bool non_uniform_n,
	const int& num_finest_elems,
	int verbose
)
{

	int xchk, ychk;
	FILE* check_fp;
	//	ChannelSegmentType* csp;  //local pointer to channel segment
	//	int chseg;	// channel segment loop counter
	int checkv;
	int version, LFv;

	NUMERIC_TYPE* dummy;

	LFv = (int(LF_VersionMajor) * 100) + (int(LF_VersionMinor) * 10) + int(LF_VersionInc);


	// #### file opening sequence before reading data

	// try and open the default checkpoint file, because if it exists it means a run was started and was halted partway.
	if ((check_fp = fopen(filenames.checkpointfilename, "rb")) == NULL) // open file and check it exists
	{
		// message so user knows file does not exist
		if (verbose == ON) printf("\nUnable to find default checkpoint file: %s\n", filenames.checkpointfilename);

		// check if user wants to load a different starting checkpoint file 
		if (states.checkfile == ON)
		{
			// try and open the file
			if ((check_fp = fopen(filenames.loadCheckpointFilename, "rb")) == NULL) // open file and check it exists
			{
				// message so user knows file does not exist
				if (verbose == ON) printf("\nUnable to find alternative starting checkpoint file: %s\n", filenames.loadCheckpointFilename);
			}
			else
			{
				// message so user knows alternative file opened ok
				if (verbose == ON) printf("\nAlternative checkpoint file opened: %s", filenames.loadCheckpointFilename);
			}
		}
	}
	else
	{
		if (states.checkfile == ON) // let user know program is using newer default checkpoint rather then alternative file
		{
			printf("\nAlternative checkpoint file: %s not opened as newer default file exists", filenames.loadCheckpointFilename);
			printf("\n - please delete: %s if you want to start from the Alternative file\n", filenames.checkpointfilename);
		}
		// message so user knows default file opened ok
		if (verbose == ON) printf("\nDefault checkpoint file opened: %s", filenames.checkpointfilename);
	}


	NUMERIC_TYPE* h = new NUMERIC_TYPE[d_assem_sol.length]();
	NUMERIC_TYPE* q = new NUMERIC_TYPE[d_non_uniform_nghbrs.length]();

	NUMERIC_TYPE* maxH = new NUMERIC_TYPE[num_finest_elems]();
	NUMERIC_TYPE* totalHtm = new NUMERIC_TYPE[num_finest_elems]();
	NUMERIC_TYPE* maxHtm = new NUMERIC_TYPE[num_finest_elems]();
	NUMERIC_TYPE* initHtm = new NUMERIC_TYPE[num_finest_elems]();

	size_t bytes_sol     = d_assem_sol.length * sizeof(NUMERIC_TYPE);
	size_t bytes_dis     = d_non_uniform_nghbrs.length * sizeof(NUMERIC_TYPE);
	size_t bytes         = num_finest_elems * sizeof(NUMERIC_TYPE);


	// #### reading data from a file if opened correctly sequence

	if (check_fp == NULL) // no file opened in initial sequence
	{
		// message so user knows no files opened and program will start from scratch
		if (verbose == ON) printf("\nNo checkpoint file opened: Starting from scratch\n");
	}
	else		 // go ahead and read data if a file was successfully opened
	{
		// message so user knows file exists and is being read in
		if (verbose == ON) printf("\nReading checkpoint file\n");

		// read checkpointing file version - for future compatibility
		fread(&checkv, sizeof(int), 1, check_fp);
		if (feof(check_fp)) {
			if (verbose == ON) printf("\nUnable to read from checkpoint file (%s) zero size: Starting from scratch\n", filenames.checkpointfilename);
		}
		else if (checkv != int(LF_CheckVersion)) {
			if (verbose == ON) printf("\nFile checkpoint version differs from current version: Starting from scratch\n");
		}
		else {
			// file exists and is not of zero length so go ahead and read data

			// check LISFLOOD-FP version
			fread(&version, sizeof(int), 1, check_fp);
			if (version != LFv) {
				if (verbose == ON) printf("\nWARNING: LISFLOOD-FP version differs from current version\n");
			}
			fread(&solver.itrn_time, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&solver.t, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&solver.itCount, sizeof(long), 1, check_fp);

			fread(&pars.MassTotal, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&solver.InitTstep, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&solver.MinTstep, sizeof(NUMERIC_TYPE), 1, check_fp);

			fread(&pars.SaveNo, sizeof(int), 1, check_fp);
			fread(&pars.SaveTotal, sizeof(NUMERIC_TYPE), 1, check_fp);

			//Check domain dimensions are as expected to prevent memory overflows/other odd errors
			fread(&xchk, sizeof(int), 1, check_fp);
			fread(&ychk, sizeof(int), 1, check_fp);
			if (xchk != pars.xsz || ychk != pars.ysz)
			{
				if (verbose == ON)
				{
					printf("\nxchk %i and xsz %i\n", xchk, pars.xsz);
					printf("ychk %i and ysz %i\n", ychk, pars.ysz);
					printf("Domain dimensions do not match those in the checkpoint file!\n");
				}
				exit(0);
			}


			fread(h, sizeof(NUMERIC_TYPE), d_assem_sol.length, check_fp);
			fread(q, sizeof(NUMERIC_TYPE), d_non_uniform_nghbrs.length, check_fp);

			if (non_uniform_n) {
				NUMERIC_TYPE* n0 = new NUMERIC_TYPE[d_assem_sol.length]();
				fread(n0, sizeof(NUMERIC_TYPE), d_assem_sol.length, check_fp);
				copy_cuda(d_assem_sol.n0, n0, bytes_sol);
			}


			fread(maxH, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
			fread(maxHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
			fread(initHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
			if (checkv > 1) fread(totalHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);

			//fread(&ChannelSegments->N_Channel_Segments,sizeof(int),1,check_fp);
			//fread(ChannelSegments,sizeof(ChannelSegmentType),ChannelSegments->N_Channel_Segments,check_fp);

//			fread(ChanMask, sizeof(int), pars.xsz * pars.ysz, check_fp);
//			fread(Arrptr->SegMask, sizeof(int), Parptr->xsz * Parptr->ysz, check_fp);

//			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
			//			fread(&Stats_Entry.mass.out, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
			//			fread(&Stats_Entry., sizeof(NUMERIC_TYPE), 1, check_fp); // QChanOut -> channel flow routine: no need
			//			fread(&Stats_Entry., sizeof(NUMERIC_TYPE), 1, check_fp); // Hds -> Water depth at the downstream exit of the model domain in meters
			//			fread(&Stats_Entry., sizeof(NUMERIC_TYPE), 1, check_fp);
			//			fread(&Stats_Entry., sizeof(NUMERIC_TYPE), 1, check_fp);

			fread(&pars.dx, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&pars.dy, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&pars.dA, sizeof(NUMERIC_TYPE), 1, check_fp);

			fread(&pars.tlx, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&pars.tly, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&pars.blx, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&pars.bly, sizeof(NUMERIC_TYPE), 1, check_fp);

			//			fread(&Stats_Entry.area, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
			//			fread(&Stats_Entry.volume, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&Stats_Collector.previous_volume, sizeof(NUMERIC_TYPE), 1, check_fp);

			//			fread(&Stats_Entry.discharge_error, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);

			//			if (checkv > 2) {
			//				fread(dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
			//				fread(dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
			//			}

			if (verbose == ON) printf("\n - Computation time so far: %.4" NUM_FMT" hrs, Sim time: %.3" NUM_FMT" of %.3" NUM_FMT" secs\n", (solver.itrn_time / C(3600.0)), solver.t, solver.Sim_Time);
			fclose(check_fp);
		}
	}


	copy_cuda(d_assem_sol.h, h, bytes_sol);
	copy_cuda(d_non_uniform_nghbrs.q, q, bytes_dis);

	copy_cuda(d_maxH, maxH, bytes);
	copy_cuda(d_totalHtm, totalHtm, bytes);
	copy_cuda(d_maxHtm, maxHtm, bytes);
	copy_cuda(d_initHtm, initHtm, bytes);

	delete[] h;
	delete[] q;

	delete[] maxH;
	delete[] totalHtm;
	delete[] maxHtm;
	delete[] initHtm;

	return;

}
