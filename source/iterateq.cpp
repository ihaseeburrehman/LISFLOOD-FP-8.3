/*
*****************************************************************************
ITERATEQ and UPDATEH
---------------------


*****************************************************************************
*/

#include "lisflood.h"
#include "rain/rain.h"
#include "sgc.h"

//-----------------------------------------------------------------------------
// ITERATE THROUGH TIME STEPS
void IterateQ(Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, SGCprams *SGCptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, const int verbose)
{
	int i, j, ptr;
	int tstep_counter = -1;   // start at -1 so that in first run through we calculate river
	int steadyCount = 0;
	NUMERIC_TYPE tstep_channel = 0; // channel timestep
	NUMERIC_TYPE Previous_t;      // previous time channel was calculated
	NUMERIC_TYPE FloodArea, dA;
	NUMERIC_TYPE discharge = C(0.0); // value of discharge for virtual gauge output
	NUMERIC_TYPE loss; //temp variable to keep track of losses since last mass interval
	NUMERIC_TYPE Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;
    DynamicRain<> dynamic_rain(Fnameptr->dynamicrainfilename, verbose);

	Solverptr->vol1 = DomainVol(Statesptr, Parptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);

	if (verbose == ON)
	{
		printf("\nStarting time steps: ");
		fflush(stdout);
	}
	Solverptr->itrn_time_now = Solverptr->itrn_time;

	// Populating Tstep variables prior to start of simulation
	if (Statesptr->adaptive_ts == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("adaptive mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->acceleration == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("acceleration mode\n\n");
		fflush(stdout);
	}
	else if (Statesptr->Roe == ON)
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("Roe mode\n\n");
		fflush(stdout);
	}
	else
	{
		if (Solverptr->t == 0)
		{
			Solverptr->Tstep = Solverptr->InitTstep;
			Solverptr->MinTstep = Solverptr->InitTstep;
		}
		if (verbose == ON) printf("non-adaptive mode\n\n");
		fflush(stdout);
	}
	if (Statesptr->SGC == ON)
	{
		// because the SGC model calculates the time step in UpdateH rather than during calcFPflow it needs to initalise 
		// SGCtmpTstep, which would usually be calculated in UpdateH
		Solverptr->Tstep = Solverptr->InitTstep;
		CalcT(Parptr, Solverptr, Arrptr);
		Solverptr->SGCtmpTstep = Solverptr->Tstep;
		if (verbose == ON) printf("SGC mode\n\n");
		fflush(stdout);
	}

	// set previous time to one timestep backwards so that first river calcs uses Solverptr->Tstep for 1st iteration
	// this is because of the way the timestep is calculated as the time difference from the last time the river was run
	//Previous_t = -Solverptr->Tstep;
	Previous_t = Solverptr->t - Solverptr->Tstep;

	time_t loop_start;
	time(&loop_start);

	// main iteration loop
	while (Solverptr->t < Solverptr->Sim_Time)
	{

		// ChannelQ routine is run using Tstep calculated for previous iteration
		// Use default kinematic solver unless diffusive flag set
		if (Statesptr->ChannelPresent == ON)
		{
			if (tstep_counter == -1 || tstep_counter == Solverptr->ts_multiple)
				// do river calc if start of run (-1) or if timestep counter is equal to ts_multiple
			{
				tstep_channel = Solverptr->t - Previous_t; // calculate river timestep (deltaT since last river calc)
				if (Statesptr->diffusive == ON) ChannelQ_Diff(tstep_channel, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, RiversIndexVecPtr, RiversIndexPtr);
				else ChannelQ(tstep_channel, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, RiversIndexVecPtr, RiversIndexPtr);
				tstep_counter = 0;  // set timestep counter to zero
				Previous_t = Solverptr->t; // record previous timestep
			}
			tstep_counter++; // increment timestep counter
		}
		// Chose between a sub-grid or conventional floodplain model
		if (Statesptr->SGC == ON)
		{
			printf("SGC not processed here\n");
			//// SGC uses time step calculated in update H and stored in Solverptr->SGCtmpTstep 
			//Solverptr->Tstep = Solverptr->SGCtmpTstep;
			//// sub grid floodplain models
			//SGC_FloodplainQ(Statesptr, Parptr, Solverptr, Arrptr, SGCptr);

			//SGC_BCs(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);
			//// Infiltration, evaporation and rainfall routines after time step update (TJF)11
			//if (Statesptr->calc_evap == ON) SGC_Evaporation(Parptr, Solverptr, Arrptr, SGCptr);
			//if (Statesptr->rainfall == ON && Statesptr->routing == OFF) SGC_Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
			//if (Statesptr->routing == ON) SGC_Routing(Statesptr, Parptr, Solverptr, Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
			//if (Statesptr->hazard == ON) UpdateV(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);
			//NUMERIC_TYPE Hmax = SGC_UpdateH(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, SGCptr);

			//// now calulate the next time step to be used if it is smaller than inital time step
			////Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->dx/(sqrt(Solverptr->g*Hmax)),Solverptr->InitTstep); //CCS_deletion
			//Solverptr->SGCtmpTstep = getmin(Solverptr->cfl*Parptr->min_dx_dy*Parptr->SGC_m / (sqrt(Solverptr->g*Hmax)), Solverptr->InitTstep);

			//BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);
			//// NOTES: Update Q's handeled within SGC flux equations (SGC_FloodplainQ etc.) time-step calculation intergrated with UpdateH
		}
		else // conventional floodplain models
		{
			// Time step is reset to initial time step at the start of each FloodplainQ calculation
			if (Solverptr->t>C(0.0)) Solverptr->Tstep = Solverptr->InitTstep;

			FloodplainQ(Statesptr, Parptr, Solverptr, Arrptr, SGCptr);

			BCs(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);
            drain_nodata_water(Parptr, Solverptr, BCptr, Arrptr);
			if (Statesptr->drychecking == ON) DryCheck(Parptr, Solverptr, Arrptr);

			// Infiltration, evaporation and rainfall routines after time step update (TJF)
			if (Statesptr->calc_infiltration == ON) FPInfiltration(Parptr, Solverptr, Arrptr);
			if (Statesptr->calc_evap == ON) Evaporation(Parptr, Solverptr, Arrptr);
			if (Statesptr->rainfall == ON && Statesptr->routing == OFF) Rainfall(Parptr, Solverptr, Arrptr); // CCS rainfall with routing scheme disabled
			if (Statesptr->routing == ON) Routing(Statesptr, Parptr, Solverptr, Arrptr);	// CCS Routing scheme (controls rainfall if enabled; can be used without rainfall)
            dynamic_rain.update_H(Parptr, Solverptr, Arrptr);

			// If Roe == ON used UpdateHRoe else use the normal UpdateH (JN/IV)
			// if(Statesptr->Roe==ON) UpdateHRoe(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
			// else UpdateH(Statesptr,Parptr,Solverptr,BCptr,ChannelSegments,Arrptr);
			// change to CalcFPQxRoe and CalcFPQyRoe... now return Q's for UpdateH (JCN)
			UpdateH(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);

			if (Statesptr->hazard == ON) UpdateV(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr);

			BoundaryFlux(Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr);
			if (Statesptr->acceleration == ON) UpdateQs(Parptr, Arrptr); // Old value of Q updated with current timestep values
			if (Statesptr->Roe == ON) UpdateQsRoe(Parptr, Solverptr, Arrptr); // If Roe==on UpdateQsRoe(JN/IV)
		}

		// Update t with final Tstep calculated
		if (Solverptr->t>C(0.0)) Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
		Solverptr->t += Solverptr->Tstep;
		Solverptr->itCount += 1;
		// Update time of initial flood inundation (in hours) and total inundation time (in seconds)
		if ((Statesptr->reset_timeinit == ON) & (Solverptr->t>Parptr->reset_timeinit_time))
		{ //reset the time of initial inundation if called for in parameter file
			Statesptr->reset_timeinit = OFF;

			//#pragma omp parallel for private(j,ptr) // irrelevent benefit JCN
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				Arrptr->initHtm[ptr] = (NULLVAL);
			}
			if (verbose == ON) printf("\n Time of initial inundation reset \n");
		}
		// Update maxH, maxHtm, totalHtm and initHtm at the mass interval if mint_hk is specifed in the .par file OR at every time step if not
		if (Solverptr->t >= Parptr->MassTotal || Statesptr->mint_hk == OFF)
		{
#pragma omp parallel for private (i, ptr)
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				if (Arrptr->H[ptr] > Solverptr->DepthThresh)
				{
					if (Arrptr->initHtm[ptr] == NULLVAL /*&& (Arrptr->H[ptr]>C(0.01)*/)
						Arrptr->initHtm[ptr] = Solverptr->t / C(3600.0);
					//if (Arrptr->H[ptr]>C(0.01)) 
					Arrptr->totalHtm[ptr] += Solverptr->Tstep / C(3600.0);
					// Update maximum water depths, and time of maximum (in hours) TDF updated only track maxH when h > DepthThresh
					if (Arrptr->H[ptr] > Arrptr->maxH[ptr])
					{
						Arrptr->maxH[ptr] = Arrptr->H[ptr];
						Arrptr->maxHtm[ptr] = Solverptr->t / C(3600.0);
					}
				}
			}
		}
		// Calculate mass balance error
		if (Solverptr->t >= Parptr->MassTotal)
		{
			Solverptr->vol2 = DomainVol(Statesptr, Parptr, ChannelSegments, Arrptr, ChannelSegmentsVecPtr); // CCS

			// calc losses for this mass interval
			loss = (Parptr->InfilTotalLoss - Parptr->InfilLoss) + (Parptr->EvapTotalLoss - Parptr->EvapLoss) - (Parptr->RainTotalLoss - Parptr->RainLoss);

			//Solverptr->Qerror=BCptr->Qin-BCptr->Qout-(Solverptr->vol2+loss-Solverptr->vol1)/Parptr->MassInt;
			// New version using VolInMT and VolOutMT
			// volume error
			Solverptr->Verror = BCptr->VolInMT - BCptr->VolOutMT - (Solverptr->vol2 + loss - Solverptr->vol1);

			// Q error
			Solverptr->Qerror = Solverptr->Verror / Parptr->MassInt;
			// reset to C(0.0)
			BCptr->VolInMT = C(0.0);
			BCptr->VolOutMT = C(0.0);

			// record cumulative loss for next time.
			Parptr->InfilLoss = Parptr->InfilTotalLoss;
			Parptr->EvapLoss = Parptr->EvapTotalLoss;
			Parptr->RainLoss = Parptr->RainTotalLoss;

			// Calculate flood area
			FloodArea = C(0.0);
			dA = Parptr->dA;
#pragma omp parallel for private(i,ptr) firstprivate(dA) reduction ( + : FloodArea)
			for (j = 0; j<Parptr->ysz; j++) for (i = 0; i<Parptr->xsz; i++)
			{
				ptr = i + j*Parptr->xsz;
				if (Statesptr->latlong == ON) dA = Arrptr->dA[ptr]; // if latlong is on change dA to local cell area
				if (Statesptr->porosity == ON)
				{
					if (Arrptr->H[ptr]>C(0.01)) FloodArea += dA*Arrptr->paerial[ptr]; // If porosity used, scale flooded area by porosity (TJF)
				}
				else if (Statesptr->SGC == ON)
				{
					if (Arrptr->H[ptr] - Arrptr->SGCbfH[ptr] > Solverptr->DepthThresh)
						FloodArea += dA; // If sub-grid used remove channel depth
				}
				else
				{
					if (Arrptr->H[ptr] > C(0.01)) FloodArea += dA; // standard area calculation
				}
			}
			Solverptr->FArea = FloodArea;

			// output mass variables line to file
			//fprintf(Fptr->mass_fp,"%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4e %12.4e  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4e %12.4e %12.4e\n",Solverptr->t,Solverptr->Tstep,Solverptr->MinTstep,Solverptr->itCount,Solverptr->FArea,Solverptr->vol2,BCptr->Qin,Solverptr->Hds,BCptr->Qout,Solverptr->Qerror,Solverptr->Verror,Parptr->RainTotalLoss-(Parptr->InfilTotalLoss+Parptr->EvapTotalLoss));
			fprintf(Fptr->mass_fp, "%-12.3" NUM_FMT" %-10.4" NUM_FMT" %-10.4" NUM_FMT" %-10ld %12.4" NUM_FMT" %12.4" NUM_FMT"  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT" %12.4" NUM_FMT"\n", Solverptr->t, Solverptr->Tstep, Solverptr->MinTstep, Solverptr->itCount, Solverptr->FArea, Solverptr->vol2, BCptr->Qin, Solverptr->Hds, BCptr->Qout, Solverptr->Qerror, Solverptr->Verror, Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss));
			fflush(Fptr->mass_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.

			Solverptr->vol1 = Solverptr->vol2;
			Parptr->MassTotal += Parptr->MassInt;

			//stage output
			if (Statesptr->save_stages == ON)
			{
				fprintf(Fptr->stage_fp, "%12.3" NUM_FMT"", Solverptr->t);
				for (i = 0; i<Locptr->Nstages; i++)
				{
					int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * Parptr->xsz;
					if (Locptr->stage_check[i] == 1) fprintf(Fptr->stage_fp, "%10.4" NUM_FMT"", Arrptr->H[index]);
					else fprintf(Fptr->stage_fp, "-\t");
				}
				fprintf(Fptr->stage_fp, "\n");
				fflush(Fptr->stage_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
				// added to export scalar velocity
				if (Statesptr->voutput_stage == ON)
				{
					fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", Solverptr->t);
					for (i = 0; i<Locptr->Nstages; i++)
					{
						int index = Locptr->stage_grid_x[i] + Locptr->stage_grid_y[i] * (Parptr->xsz + 1);
						if (Locptr->stage_check[i] == 1) fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", sqrt(pow(getmax(fabs(Arrptr->Vx[index]), fabs(Arrptr->Vx[index + 1])), 2) + pow(getmax(fabs(Arrptr->Vy[index]), fabs(Arrptr->Vy[index + (Parptr->xsz + 1)])), 2)));
						else fprintf(Fptr->vel_fp, "-\t");
					}
					fprintf(Fptr->vel_fp, "\n");
					fflush(Fptr->vel_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
				}
			}

			//virtual gauge output
			if (Statesptr->gsection == ON)
			{
				fprintf(Fptr->gau_fp, "%12.2" NUM_FMT"", Solverptr->t); // print tiem to file
				for (i = 0; i<Locptr->Ngauges; i++) // loop through each virtual gauge
				{
					// call discharge calculation function
					discharge = CalcVirtualGauge(i, Parptr->xsz + 1, Arrptr->Qx, Arrptr->Qy, Locptr);

					fprintf(Fptr->gau_fp, " %10.3" NUM_FMT"", discharge); // Print discharge to file
				}
				fprintf(Fptr->gau_fp, "\n"); // print end of line
				fflush(Fptr->gau_fp); // force program to flush buffer to file - keeps file in sync with writes - user sometimes tracks progress through the file.
			}

			// Checkpointing
			if (Statesptr->checkpoint == ON)
			{
				//iteration time
				time(&Solverptr->time_check);
				Solverptr->itrn_time_now = Solverptr->itrn_time + (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start);
				if (Solverptr->itrn_time_now >= Parptr->nextcheck)
				{
					WriteCheckpoint(Fnameptr, Statesptr, Parptr, Solverptr, BCptr, ChannelSegments, Arrptr, verbose);
					Parptr->nextcheck = Solverptr->itrn_time_now + (Parptr->checkfreq * 3600);
				}
			}
		}

		// Regular output

		if (Solverptr->t >= Parptr->SaveTotal)
		{

			time(&Solverptr->time_check);
			Comp_time = (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start) / 60;
			if (Comp_time != 0 && Statesptr->comp_out == ON) // only of t is not zero (can't divide by zero)
			{
				Model_Comp_Ratio = ((Solverptr->t / 60) / Comp_time);
				Model_time_left = (Solverptr->Sim_Time - Solverptr->t) / 60;
				Est_Time_Fin = (Model_time_left / Model_Comp_Ratio);
				Est_Time_Tot = Comp_time + Est_Time_Fin;
				printf("T(mins): M: %.1" NUM_FMT", C: %.1" NUM_FMT", M/C: %.2" NUM_FMT", ETot: %.1" NUM_FMT", EFin: %.1" NUM_FMT"\n", (Solverptr->t / C(60.0)), Comp_time, Model_Comp_Ratio, Est_Time_Tot, Est_Time_Fin);
			}

			write_regular_output(Fnameptr, Solverptr, Statesptr, Parptr, Arrptr, SGCptr);

			//regular profile output, if requested in param file
			if (Statesptr->profileoutput == ON)
			{
				write_profile(Fnameptr->resrootname, Parptr->SaveNo, ".profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr); // CCS
			}

			// update interval counter
			Parptr->SaveTotal += Parptr->SaveInt;
			Parptr->SaveNo += 1;

		}

		// Single overpass
		if (Statesptr->single_op == ON && Solverptr->t >= Parptr->op)
		{
			// write rasters
			if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opb", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
			else write_ascfile(Fnameptr->resrootname, -1, ".op", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
			// write profiles
			if (Statesptr->profileoutput == ON) write_profile(Fnameptr->resrootname, -1, ".profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr);
			if (verbose == ON) printf("Writing overpass at %" NUM_FMT" seconds\n", Solverptr->t);
			Statesptr->single_op = OFF;
			// raster elevation output
			if (Statesptr->save_elev == ON)
			{
				if (Statesptr->SGC == ON) // SGC output
				{
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opelevb", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
					else write_ascfile(Fnameptr->resrootname, -1, ".opelev", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
				}
				else // standard model output
				{
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".opelevb", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
					else write_ascfile(Fnameptr->resrootname, -1, ".opelev", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
				}
			}
		}

		// Multiple overpasses
		if (Statesptr->multi_op == ON)
		{
			for (i = 0; i<Parptr->op_multinum; i++)
			{
				if (Solverptr->t >= Parptr->op_multisteps[i] && Parptr->op_multiswitch[i] == 0)
				{
					// raster depth output
					if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opb", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
					else write_ascfile(Fnameptr->resrootname, i, "-T.op", Arrptr->H, Arrptr->DEM, 0, Statesptr, Parptr);
					// profiles output
					if (Statesptr->profileoutput == ON) write_profile(Fnameptr->resrootname, i, "-T.profile", Statesptr, ChannelSegments, Arrptr, Parptr, RiversIndexVecPtr, RiversIndexPtr);

					Parptr->op_multiswitch[i] = 1;
					if (verbose == ON) printf("Writing overpass %d at %" NUM_FMT" seconds\n", i, Parptr->op_multisteps[i]);
					// raster elevation output
					if (Statesptr->save_elev == ON)
					{
						if (Statesptr->SGC == ON) // SGC output
						{
							if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opelevb", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
							else write_ascfile(Fnameptr->resrootname, i, "-T.opelev", Arrptr->H, Arrptr->SGCz, 3, Statesptr, Parptr, Solverptr->DepthThresh);
						}
						else // standard model output
						{
							if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, i, "-T.opelevb", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
							else write_ascfile(Fnameptr->resrootname, i, "-T.opelev", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr, Solverptr->DepthThresh);
						}
					}
				}
			}
		}

		// If requested on command line, check whether we should kill this simulation...
		if (Statesptr->killsim == ON)
		{
			//iteration time
			time(&Solverptr->time_check);
			Solverptr->itrn_time_now = Solverptr->itrn_time + (NUMERIC_TYPE)difftime(Solverptr->time_check, Solverptr->time_start);
			// check if we have reached the kill time
			if (Solverptr->itrn_time_now >= Parptr->killsim_time) {
				if (verbose == ON) printf("Simulation kill time reached... ");
				break;
			}
		}

		// If an auto steady-state has been requested:
		if (Statesptr->steadycheck == ON && Solverptr->t >= Parptr->steadyTotal)
		{
			Parptr->steadyQdiff = BCptr->Qin - BCptr->Qout;
			if (fabs(Parptr->steadyQdiff)<Parptr->steadyQtol) {
				steadyCount += 1; // keep going a bit to make sure...
				if (steadyCount == 10) break;
			}
			else {
				steadyCount = 0;
			}
			Parptr->steadyTotal += Parptr->steadyInt;
		}

	}
	//END main ITERATIONS

	time_t loop_end;
	time(&loop_end);

	double seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);


	//output regular files
	fileoutput(Fnameptr, Statesptr, Parptr, Arrptr);

	// Must be last because maxH changed to Maximum elevation !!!
	// Write maximum elevation
	if (Statesptr->SGC == ON) // SGC output
	{
		for (i = 0; i<Parptr->xsz; i++) for (j = 0; j<Parptr->ysz; j++)
		{
			ptr = i + j*Parptr->xsz;
			if (Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr] += Arrptr->SGCz[ptr]; else Arrptr->maxH[ptr] = NULLVAL;
		}
		if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".mxeb", Arrptr->maxH, Arrptr->SGCz, 0, Statesptr, Parptr);
		else write_ascfile(Fnameptr->resrootname, -1, ".mxe", Arrptr->maxH, Arrptr->SGCz, 0, Statesptr, Parptr);
	}
	else
	{
		for (i = 0; i<Parptr->xsz; i++) for (j = 0; j<Parptr->ysz; j++)
		{
			ptr = i + j*Parptr->xsz;
			if (Arrptr->maxH[ptr]>1e-3) Arrptr->maxH[ptr] += Arrptr->DEM[ptr]; else Arrptr->maxH[ptr] = NULLVAL;
		}
		if (Statesptr->binary_out == ON) write_binrasterfile(Fnameptr->resrootname, -1, ".mxeb", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
		else write_ascfile(Fnameptr->resrootname, -1, ".mxe", Arrptr->maxH, Arrptr->DEM, 0, Statesptr, Parptr);
	}

	if (verbose == ON) printf("Finished.\n\n");

	return;
}

//----------------------------------------------------------------------------
// SUM Qs INTO A CELL AND UPDATE DEPTH ACCORDINGLY
void UpdateH(States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr)
{
	int i, j;
	NUMERIC_TYPE *qxptr0, *qyptr0, *qyptr1, *hptr;
	int *mptr;
	NUMERIC_TYPE dV, himp, qtmp;
	NUMERIC_TYPE dAPorTemp;
	NUMERIC_TYPE Qpnt;

	// Insert point sources ((MT) moved before dV check, in case flow out of cell results in negative H reset and the inflow would prevent this)
	// H point sources moved back to after UpdateH (JCN)
	BCptr->Qpoint_pos = C(0.0);
	BCptr->Qpoint_neg = C(0.0);

	for (int ps_index = 0; ps_index<BCptr->numPS; ps_index++)
	{
		if (BCptr->PS_Ident[ps_index] == QFIX4) // QFIX
		{
			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] += BCptr->PS_Val[ps_index] * Parptr->dx*Solverptr->Tstep / Parptr->dA;
			Qpnt = BCptr->PS_Val[ps_index] * Parptr->dx;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;
		}
		if (BCptr->PS_Ident[ps_index] == QVAR5) // QVAR
		{
			qtmp = InterpolateTimeSeries(BCptr->PS_TimeSeries[ps_index], Solverptr->t);
			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] += qtmp*Parptr->dx*Solverptr->Tstep / Parptr->dA;
			Qpnt = qtmp*Parptr->dx;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;
		}
	}

	// Calculate dV (-ve => outflow) and update NewH
#pragma omp parallel for private( i,qxptr0,qyptr0,qyptr1,hptr,mptr,dV,dAPorTemp)
	for (j = 0; j<Parptr->ysz; j++)
	{
		qxptr0 = Arrptr->Qx + j*(Parptr->xsz + 1);
		qyptr0 = Arrptr->Qy + j*(Parptr->xsz + 1);
		qyptr1 = Arrptr->Qy + (j + 1)*(Parptr->xsz + 1);
		hptr = Arrptr->H + j*Parptr->xsz;
		mptr = Arrptr->ChanMask + j*Parptr->xsz;
		for (i = 0; i<Parptr->xsz; i++)
		{
			if (*mptr == -1)
			{
				if (Statesptr->porosity == ON)
				{
					dV = Solverptr->Tstep*(*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
					dAPorTemp = PorArea(i, j, Parptr, Arrptr);
					if (dAPorTemp == C(0.0)) (*hptr) += C(0.0);
					else (*hptr) += dV / dAPorTemp;

					if (*hptr<C(0.0)) *hptr = C(0.0);
				}
				else
				{
					dV = Solverptr->Tstep*(*qxptr0 - *(qxptr0 + 1) + *qyptr0 - *qyptr1);
					(*hptr) += dV / Parptr->dA;
					if (*hptr<C(0.0)) *hptr = C(0.0);


					//if (j == 0 && i == 79) {

					//	bool test;
					//	test = 1;
					//}
				}
			}
			qxptr0++;
			qyptr0++;
			qyptr1++;
			hptr++;
			mptr++;
		}
	}

	// Point source HVAR and HFIX
	for (int ps_index = 0; ps_index<BCptr->numPS; ps_index++)
	{
		if (BCptr->PS_Ident[ps_index] == HFIX2) // HFIX
		{
			himp = BCptr->PS_Val[ps_index] - Arrptr->DEM[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz];
			if (himp<C(0.0)) himp = C(0.0);

			Qpnt = (himp - Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz])*Parptr->dA / Solverptr->Tstep;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;

			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] = himp;
		}
		if (BCptr->PS_Ident[ps_index] == HVAR3) // HVAR
		{
			himp = InterpolateTimeSeries(BCptr->PS_TimeSeries[ps_index], Solverptr->t) - Arrptr->DEM[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz];
			if (himp<C(0.0)) himp = C(0.0);

			Qpnt = (himp - Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz])*Parptr->dA / Solverptr->Tstep;
			if (Qpnt>0) BCptr->Qpoint_pos += Qpnt;
			else BCptr->Qpoint_neg += Qpnt;

			Arrptr->H[BCptr->xpi[ps_index] + BCptr->ypi[ps_index] * Parptr->xsz] = himp;
		}
	}

	return;
}
