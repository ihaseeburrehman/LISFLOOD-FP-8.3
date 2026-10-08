#include "output.h"

void write_mass_stats
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	BoundCs *BCptr
)
{
	fprintf
	(
		Fptr->mass_fp,
		"%-12.3f %10.4e %10.4e %-10li %12.4e %12.4e  %-11.3" NUM_FMT" %-10.3" NUM_FMT" %-11.3" NUM_FMT" %12.4e %12.4e %12.4e\n",
		Solverptr->t,
		Solverptr->Tstep,
		Solverptr->MinTstep,
		Solverptr->itCount,
		Solverptr->FArea,
		Solverptr->vol2,
		BCptr->Qin,
		Solverptr->Hds,
		BCptr->Qout,
		Solverptr->Qerror,
		Solverptr->Verror,
		Parptr->RainTotalLoss - (Parptr->InfilTotalLoss + Parptr->EvapTotalLoss)
	);
	fflush(Fptr->mass_fp);
}

void write_depth_samples
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	Stage *Stageptr,
	Arrays *Arrptr
)
{
	fprintf(Fptr->stage_fp, "%12.3" NUM_FMT"", Solverptr->t);
	for (int i = 0; i<Stageptr->Nstages; i++)
	{
		NUMERIC_TYPE H = Arrptr->H[
			Stageptr->stage_grid_y[i]*Parptr->xsz + Stageptr->stage_grid_x[i]];

		if (Stageptr->stage_check[i] == 1)
		{
			fprintf(Fptr->stage_fp, "%10.4" NUM_FMT"", H);
		}
		else
		{
			fprintf(Fptr->stage_fp, "-\t");
		}
	}
	fprintf(Fptr->stage_fp, "\n");
	fflush(Fptr->stage_fp);
}

void write_speed_samples
(
	Files *Fptr,
	Pars *Parptr,
	Solver *Solverptr,
	Stage *Stageptr,
	Arrays *Arrptr
)
{
	fprintf(Fptr->vel_fp, "%12.3" NUM_FMT"", Solverptr->t);
	for (int i = 0; i<Stageptr->Nstages; i++)
	{
		const int index = 
			Stageptr->stage_grid_y[i]*Parptr->xsz + Stageptr->stage_grid_x[i];

		NUMERIC_TYPE H = Arrptr->H[index];
		NUMERIC_TYPE HU = Arrptr->HU[index];
		NUMERIC_TYPE HV = Arrptr->HV[index];

		NUMERIC_TYPE speed;

		if (H > Solverptr->DepthThresh)
		{
			speed = sqrt(HU/H * HU/H + HV/H * HV/H);
		}
		else
		{
			speed = C(0.0);
		}

		if (Stageptr->stage_check[i] == 1)
		{
			fprintf(Fptr->vel_fp, "%10.4" NUM_FMT"", speed);
		}
		else
		{
			fprintf(Fptr->vel_fp, "-\t");
		}
	}
	fprintf(Fptr->vel_fp, "\n");
	fflush(Fptr->vel_fp);
}

void write_solution
(
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr,
	Solver *Solverptr,
	Arrays *Arrptr
)
{
	if (Solverptr->t < Parptr->SaveTotal) return;

    if (Statesptr->save_depth == ON)
	{
		write_field(Arrptr->H, ".wd", ".wdb", Fnameptr, Statesptr, Parptr);
	}

	if (Statesptr->save_elev == ON)
	{
		if (Statesptr->binary_out == ON)
		{
			write_binrasterfile(Fnameptr->resrootname, Parptr->SaveNo,
					".elevb", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr);
		}
		else
		{
			write_ascfile(Fnameptr->resrootname, Parptr->SaveNo,
					".elev", Arrptr->H, Arrptr->DEM, 3, Statesptr, Parptr);
		}
	}
	
    if (Statesptr->save_Qs == ON)
	{
		write_field(Arrptr->HU, ".Qx", ".Qxb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->HV, ".Qy", ".Qyb", Fnameptr, Statesptr, Parptr);
	}

	if (Statesptr->voutput == ON)
	{
		write_field(Arrptr->Vx, ".Vx", ".Vxb", Fnameptr, Statesptr, Parptr);
		write_field(Arrptr->Vy, ".Vy", ".Vyb", Fnameptr, Statesptr, Parptr);
	}

	Parptr->SaveTotal += Parptr->SaveInt;
	Parptr->SaveNo += 1;
}

void write_max_field
(
	Fnames* Fnameptr,
	States* Statesptr,
    Pars* Parptr,
    Arrays* Arrptr
)
{
	if (Statesptr->binary_out == ON)
	{
	    write_binrasterfile(Fnameptr->resrootname, -1, ".maxb", Arrptr->maxH,
                Arrptr->DEM, 0, Statesptr, Parptr);
    }
    else
    {
	    write_ascfile(Fnameptr->resrootname, -1, ".max", Arrptr->maxH,
                Arrptr->DEM, 0, Statesptr, Parptr);
    }
}

void write_field
(
	NUMERIC_TYPE *field,
	const char *ascii_suffix,
	const char *binary_suffix,
	Fnames *Fnameptr,
	States *Statesptr,
	Pars *Parptr
)
{
	if (Statesptr->binary_out == ON)
	{
		write_binrasterfile(Fnameptr->resrootname, Parptr->SaveNo,
				binary_suffix, field, nullptr, 0, Statesptr, Parptr);
	}
	else
	{
		char format_specifier[16];
		snprintf(format_specifier, 16*sizeof(char), "%%.%i",
				Parptr->output_precision);

		write_ascfile(Fnameptr->resrootname, Parptr->SaveNo,
				ascii_suffix, field, nullptr, 0, Statesptr, Parptr, C(0.0), format_specifier);
	}
}

