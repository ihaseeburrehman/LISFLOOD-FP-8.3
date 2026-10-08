#include "write_chkpnt.cuh"
#include "../../VersionHistory.h"

__host__ void lis::cuda::acc_nugrid::write_chkpnt
(
	const AssembledSolution&    d_assem_sol,
	const NonUniformNeighbours& d_non_uniform_nghbrs,
	const Fnames& filenames,
	States& states,
	const Pars& pars,
	const Solver& solver,
	NUMERIC_TYPE* d_maxH,
	NUMERIC_TYPE* d_totalHtm,
	NUMERIC_TYPE* d_maxHtm,
	NUMERIC_TYPE* d_initHtm,
	bool non_uniform_n,
	const int& num_finest_elems,
	int verbose
)
{

	FILE* check_fp;
	int checkv;
	int version;

	checkv = int(LF_CheckVersion);
	version = (int(LF_VersionMajor) * 100) + (int(LF_VersionMinor) * 10) + int(LF_VersionInc);

	NUMERIC_TYPE* h = new NUMERIC_TYPE[d_assem_sol.length];
	NUMERIC_TYPE* q = new NUMERIC_TYPE[d_non_uniform_nghbrs.length];

	size_t bytes_sol     = d_assem_sol.length * sizeof(NUMERIC_TYPE);
	size_t bytes_dis     = d_non_uniform_nghbrs.length * sizeof(NUMERIC_TYPE);

	copy_cuda(h, d_assem_sol.h, bytes_sol);
	copy_cuda(q, d_non_uniform_nghbrs.q, bytes_dis);

	NUMERIC_TYPE* maxH = new NUMERIC_TYPE[num_finest_elems];
	NUMERIC_TYPE* totalHtm = new NUMERIC_TYPE[num_finest_elems];
	NUMERIC_TYPE* maxHtm = new NUMERIC_TYPE[num_finest_elems];
	NUMERIC_TYPE* initHtm = new NUMERIC_TYPE[num_finest_elems];

	size_t bytes = num_finest_elems * sizeof(NUMERIC_TYPE);

	copy_cuda(maxH, d_maxH, bytes);
	copy_cuda(totalHtm, d_totalHtm, bytes);
	copy_cuda(maxHtm, d_maxHtm, bytes);
	copy_cuda(initHtm, d_initHtm, bytes);

	//File written in binary rather than ASCII to maintain model precision
	if ((check_fp = fopen(filenames.checkpointfilename, "wb")) == NULL) {
		if (verbose == ON) printf("Unable to open checkpoint file: checkpointing off");
		states.checkpoint = OFF;
		return;
	}
	else {
		//checkpointing file version - for future compatibility
		fwrite(&checkv, sizeof(int), 1, check_fp);
		fwrite(&version, sizeof(int), 1, check_fp);

		fwrite(&solver.itrn_time_now, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&solver.t, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&solver.itCount, sizeof(long), 1, check_fp);

		fwrite(&pars.MassTotal, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&solver.Tstep, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&solver.MinTstep, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&pars.SaveNo, sizeof(int), 1, check_fp);
		fwrite(&pars.SaveTotal, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&pars.xsz, sizeof(int), 1, check_fp);
		fwrite(&pars.ysz, sizeof(int), 1, check_fp);

		fwrite(h, sizeof(NUMERIC_TYPE), d_assem_sol.length, check_fp);
		fwrite(q, sizeof(NUMERIC_TYPE), d_non_uniform_nghbrs.length, check_fp);

		if (non_uniform_n) {
			NUMERIC_TYPE* n0 = new NUMERIC_TYPE[d_assem_sol.length];
			copy_cuda(n0, d_assem_sol.n0, bytes_sol);
			fwrite(n0, sizeof(NUMERIC_TYPE), d_assem_sol.length, check_fp);
		}

		fwrite(maxH, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
		fwrite(maxHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
		fwrite(initHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);
		if (checkv > 1) fwrite(totalHtm, sizeof(NUMERIC_TYPE), num_finest_elems, check_fp);

//		fwrite(&Stats_Entry.mass.in, sizeof(NUMERIC_TYPE), 1, check_fp);
//		fwrite(&Stats_Entry.mass.out, sizeof(NUMERIC_TYPE), 1, check_fp);
		//		fwrite(&BCptr->QChanOut, sizeof(NUMERIC_TYPE), 1, check_fp);
		//		fwrite(&Solverptr->Hds, sizeof(NUMERIC_TYPE), 1, check_fp);
		//		fwrite(&BCptr->Qpoint_pos, sizeof(NUMERIC_TYPE), 1, check_fp);
		//		fwrite(&BCptr->Qpoint_neg, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&pars.dx, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&pars.dy, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&pars.dA, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&pars.tlx, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&pars.tly, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&pars.blx, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&pars.bly, sizeof(NUMERIC_TYPE), 1, check_fp);

//		fwrite(&Stats_Entry.area, sizeof(NUMERIC_TYPE), 1, check_fp);
		//		fwrite(&Solverptr->vol1, sizeof(NUMERIC_TYPE), 1, check_fp);
//		fwrite(&Stats_Entry.volume, sizeof(NUMERIC_TYPE), 1, check_fp);
//		fwrite(&Stats_Entry.discharge_error, sizeof(NUMERIC_TYPE), 1, check_fp);

		//if (checkv > 2) {
		//	fwrite(&Parptr->InfilTotalLoss, sizeof(NUMERIC_TYPE), 1, check_fp);
		//	fwrite(&Parptr->EvapTotalLoss, sizeof(NUMERIC_TYPE), 1, check_fp);
		//}

		fclose(check_fp);
		printf("Checkpointed at %.4" NUM_FMT" hours computation time\n", (solver.itrn_time_now / C(3600.0)));
	}


	delete[] h;
	delete[] q;

	delete[] maxH;
	delete[] totalHtm;
	delete[] maxHtm;
	delete[] initHtm;

	return;

}
