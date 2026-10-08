#include "cuda_simulate.cuh"
#include "cuda_util.cuh"
#include "cuda_unifiedallocator.cuh"
#include "io.h"
#include "../VersionHistory.h"
//#include "cuda_stats.cuh"
//#include "acc/cuda_acc_flow.cuh"

void lis::cuda::Simulation::print_device_info()
{
	int device = cuda::get_device();
	cudaDeviceProp cudaDevProp;
	cuda::get_device_properties(cudaDevProp, device);
	printf("Using CUDA device: %i\t", device);
	unsigned char* bytes = reinterpret_cast<unsigned char*>(cudaDevProp.uuid.bytes);
	printf("UUID: GPU-%02x%02x%02x%02x-", bytes[0], bytes[1], bytes[2], bytes[3]);
	printf("%02x%02x-", bytes[4], bytes[5]);
	printf("%02x%02x-", bytes[6], bytes[7]);
	printf("%02x%02x-", bytes[8], bytes[9]);
	printf("%02x%02x%02x%02x%02x%02x\n", bytes[10], bytes[11], bytes[12], bytes[13], bytes[14], bytes[15]);
}

void lis::cuda::Simulation::initialise_H
(
	NUMERIC_TYPE* H,
	const char* filename,
	States& states,
	NUMERIC_TYPE* DEM,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	if (states.startfile == ON)
	{
		StartFile::load(filename, H, geometry, pitch, offset, verbose);
		if (states.startelev == ON)
		{
			StartFile::subtract_dem(H, DEM, geometry, pitch, offset);
		}
	}


}

void lis::cuda::Simulation::nullify_max
(
	NUMERIC_TYPE* maxHtm,
	NUMERIC_TYPE* initHtm,
	NUMERIC_TYPE* totalHtm,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{

	for (int j = 0; j < geometry.ysz; j++)
	{
		for (int i = 0; i < geometry.xsz; i++)
		{
			maxHtm[j * pitch + i + offset] = NULLVAL;
			initHtm[j * pitch + i + offset] = NULLVAL;
			totalHtm[j * pitch + i + offset] = C(0.0);

		}
	}
}

void lis::cuda::Simulation::initialise_discharge
(
	NUMERIC_TYPE* HU,
	NUMERIC_TYPE* HV,
	const char* startfilename,
	States& states,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	if (states.startq2d == ON)
	{
		char hu_startfile[800];
		strcpy(hu_startfile, startfilename);
		strcat(hu_startfile, ".Qx");

		char hv_startfile[800];
		strcpy(hv_startfile, startfilename);
		strcat(hv_startfile, ".Qy");

		StartFile::load(hu_startfile, HU, geometry, pitch, offset, verbose);
		StartFile::load(hv_startfile, HV, geometry, pitch, offset, verbose);
	}
}

void lis::cuda::Simulation::initialise_manning
(
	NUMERIC_TYPE* manning,
	const char* filename,
	Geometry& geometry,
	int pitch,
	int offset,
	int verbose
)
{
	FILE* file = fopen_or_die(filename, "rb", "Loading manningfile\n", verbose);
	Geometry manning_file_geometry;
	NUMERIC_TYPE no_data_value;
	AsciiRaster::read_header(file, manning_file_geometry, no_data_value);
	AsciiRaster::match_cell_dimensions_or_die(geometry, manning_file_geometry,
			"initialise_manning");
	AsciiRaster::read(file, manning, geometry, pitch, offset);
	fclose(file);
}

void lis::cuda::Simulation::load_boundaries
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	BoundCs& boundCs,
	int verbose
)
{
	LoadBCs(&filenames, &states, &pars, &boundCs, verbose);
	LoadBCVar(&filenames, &states, &pars, &boundCs, nullptr, nullptr, nullptr,
			verbose);
}

void lis::cuda::Simulation::update_geometry
(
	Pars& dst,
	Geometry& src
)
{
	dst.xsz = src.xsz;
	dst.ysz = src.ysz;
	dst.blx = src.blx;
	dst.bly = src.bly;
	dst.tly = src.tly;
	dst.dx = src.dx;
	dst.dy = src.dy;
}

//read_checkpoint(filenames, states, pars, solver, BCptr, CSTypePtr, Arrptr, verbose);

void lis::cuda::Simulation::read_checkpoint
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	Solver& solver,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* Qx,
	NUMERIC_TYPE* Qy,
	NUMERIC_TYPE* maxH,
	NUMERIC_TYPE* totalHtm,
	NUMERIC_TYPE* maxHtm,
	NUMERIC_TYPE* initHtm,
	StatsCollector& Stats_Collector,
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

			fread(&states.single_op, sizeof(int), 1, check_fp);
			fread(&pars.op_multiswitch, sizeof(int), pars.op_multinum, check_fp);
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

			//Malloc dependent variables - see LoadDEM
			fread(H, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
			fread(Qx, sizeof(NUMERIC_TYPE), (pars.xsz + 1) * (pars.ysz + 1), check_fp);
			fread(Qy, sizeof(NUMERIC_TYPE), (pars.xsz + 1) * (pars.ysz + 1), check_fp);
			fread(maxH, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
			fread(maxHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
			fread(initHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
			if (checkv > 1) fread(totalHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);

			//fread(&ChannelSegments->N_Channel_Segments,sizeof(int),1,check_fp);
			//fread(ChannelSegments,sizeof(ChannelSegmentType),ChannelSegments->N_Channel_Segments,check_fp);

//			fread(ChanMask, sizeof(int), pars.xsz * pars.ysz, check_fp);
//			fread(Arrptr->SegMask, sizeof(int), Parptr->xsz * Parptr->ysz, check_fp);

			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&Stats_Entry.mass.out, sizeof(NUMERIC_TYPE), 1, check_fp);
 			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
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
			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
//			fread(&Stats_Entry.volume, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&Stats_Collector.previous_volume, sizeof(NUMERIC_TYPE), 1, check_fp);

//			fread(&Stats_Entry.discharge_error, sizeof(NUMERIC_TYPE), 1, check_fp);
			fread(&dummy, sizeof(NUMERIC_TYPE), 1, check_fp);

//			if (checkv > 2) {
//				fread(dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
//				fread(dummy, sizeof(NUMERIC_TYPE), 1, check_fp);
//			}

			if (verbose == ON) printf("\n - Computation time so far: %.4" NUM_FMT" hrs, Sim time: %.3" NUM_FMT" of %.3" NUM_FMT" secs\n", (solver.itrn_time / C(3600.0)), solver.t, solver.Sim_Time);
			fclose(check_fp);
		}
	}

	return;

}


void lis::cuda::Simulation::write_checkpoint
(
	Fnames& filenames,
	States& states,
	Pars& pars,
	Solver& solver,
	NUMERIC_TYPE* H,
	NUMERIC_TYPE* Qx,
	NUMERIC_TYPE* Qy,
	NUMERIC_TYPE* maxH,
	NUMERIC_TYPE* totalHtm,
	NUMERIC_TYPE* maxHtm,
	NUMERIC_TYPE* initHtm,
	lis::StatsEntry& Stats_Entry,
	NUMERIC_TYPE dt,
	int verbose
)
{

	FILE* check_fp;
//	ChannelSegmentType* csp;  //local pointer to channel segment
//	int chseg;	// channel segment loop counter
	int checkv;
	int version;

	checkv = int(LF_CheckVersion);
	version = (int(LF_VersionMajor) * 100) + (int(LF_VersionMinor) * 10) + int(LF_VersionInc);



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
		fwrite(&dt, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&solver.MinTstep, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&states.single_op, sizeof(int), 1, check_fp);
		fwrite(pars.op_multiswitch, sizeof(int), pars.op_multinum, check_fp);
		fwrite(&pars.SaveNo, sizeof(int), 1, check_fp);
		fwrite(&pars.SaveTotal, sizeof(NUMERIC_TYPE), 1, check_fp);

		fwrite(&pars.xsz, sizeof(int), 1, check_fp);
		fwrite(&pars.ysz, sizeof(int), 1, check_fp);

		fwrite(H, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
		fwrite(Qx, sizeof(NUMERIC_TYPE), (pars.xsz + 1) * (pars.ysz + 1), check_fp);
		fwrite(Qy, sizeof(NUMERIC_TYPE), (pars.xsz + 1) * (pars.ysz + 1), check_fp);
		fwrite(maxH, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
		fwrite(maxHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
		fwrite(initHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);
		if (checkv > 1) fwrite(totalHtm, sizeof(NUMERIC_TYPE), pars.xsz * pars.ysz, check_fp);

		//fwrite(&ChannelSegments->N_Channel_Segments,sizeof(int),1,check_fp);
		//fwrite(ChannelSegments,sizeof(ChannelSegmentType),ChannelSegments->N_Channel_Segments,check_fp);

//		fwrite(ChanMask, sizeof(int), pars.xsz * pars.ysz, check_fp);
//		fwrite(Arrptr->SegMask, sizeof(int), Parptr->xsz * Parptr->ysz, check_fp);

		fwrite(&Stats_Entry.mass.in, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&Stats_Entry.mass.out, sizeof(NUMERIC_TYPE), 1, check_fp);
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

		fwrite(&Stats_Entry.area, sizeof(NUMERIC_TYPE), 1, check_fp);
//		fwrite(&Solverptr->vol1, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&Stats_Entry.volume, sizeof(NUMERIC_TYPE), 1, check_fp);
		fwrite(&Stats_Entry.discharge_error, sizeof(NUMERIC_TYPE), 1, check_fp);

		//if (checkv > 2) {
		//	fwrite(&Parptr->InfilTotalLoss, sizeof(NUMERIC_TYPE), 1, check_fp);
		//	fwrite(&Parptr->EvapTotalLoss, sizeof(NUMERIC_TYPE), 1, check_fp);
		//}

		fclose(check_fp);
		printf("Checkpointed at %.4" NUM_FMT" hours computation time\n", (solver.itrn_time_now / C(3600.0)));
	}

	return;

}