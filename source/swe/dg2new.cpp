#include "dg2new.h"

#include "dg2/dg2_output.h"
#include "hll.h"
#include "input.h"
#include "output.h"
#include "stats.h"
#include "../utility.h"
#include <algorithm>

dg2new::DG2Solver::DG2Solver
(
    Fnames *Fnameptr,
    Files *Fptr,
    States *Statesptr,
    Pars *Parptr,
    Solver *Solverptr,
    BoundCs *BCptr,
    Stage *Stageptr,
    Arrays *Arrptr,
    int verbose
)
:
Fnameptr(Fnameptr),
Fptr(Fptr),
Statesptr(Statesptr),
Parptr(Parptr),
Solverptr(Solverptr),
BCptr(BCptr),
Stageptr(Stageptr),
Arrptr(Arrptr),
verbose(verbose)
{
    read_dem_slopes();
    read_initial_conditions();
    malloc_new_fields();

	for (int j=0; j<Parptr->ysz; j++)
	{
        {
            int i = 0;
            Arrptr->DEM[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1x[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1y[j*Parptr->xsz+i] = C(0.0);
            Arrptr->H[j*Parptr->xsz+i] = C(0.0);
        }

        {
            int i = Parptr->xsz-1;
            Arrptr->DEM[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1x[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1y[j*Parptr->xsz+i] = C(0.0);
            Arrptr->H[j*Parptr->xsz+i] = C(0.0);
        }
    }

	for (int i=0; i<Parptr->xsz; i++)
	{
        {
            int j = 0;
            Arrptr->DEM[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1x[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1y[j*Parptr->xsz+i] = C(0.0);
            Arrptr->H[j*Parptr->xsz+i] = C(0.0);
        }

        {
            int j = Parptr->ysz-1;
            Arrptr->DEM[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1x[j*Parptr->xsz+i] = C(0.0);
            Arrptr->DEM1y[j*Parptr->xsz+i] = C(0.0);
            Arrptr->H[j*Parptr->xsz+i] = C(0.0);
        }
    }
}

dg2new::DG2Solver::~DG2Solver()
{
	memory_free_legacy(&H_new);
	memory_free_legacy(&HU_new);
	memory_free_legacy(&HV_new);
	memory_free_legacy(&H1x_new);
	memory_free_legacy(&HU1x_new);
	memory_free_legacy(&HV1x_new);
	memory_free_legacy(&H1y_new);
	memory_free_legacy(&HU1y_new);
	memory_free_legacy(&HV1y_new);
}

void dg2new::DG2Solver::solve()
{
    set_first_Tstep();

	time_t loop_start;
	time(&loop_start);

	while (Solverptr->t < Solverptr->Sim_Time)
	{	
        print_Tstep();

        apply_friction();
        rk_stage1();
        rk_stage2();


        update_Tstep();
        update_mass_stats();

        write_solution();
    }

	time_t loop_end;
	time(&loop_end);

	double seconds = difftime(loop_end, loop_start);
	printf("loop time %lf\n", seconds);

	if (verbose == ON) printf("Finished.\n\n");
}

void dg2new::DG2Solver::read_dem_slopes()
{
	char dem1x[256];
	char dem1y[256];
	strcpy(dem1x, Fnameptr->demfilename);
	strcat(dem1x, "1x");
	strcpy(dem1y, Fnameptr->demfilename);
	strcat(dem1y, "1y");
	read_ascfile(dem1x, Parptr, Arrptr->DEM1x, "Loading DEM1x\n", verbose);
	read_ascfile(dem1y, Parptr, Arrptr->DEM1y, "Loading DEM1y\n", verbose);
}

void dg2new::DG2Solver::read_initial_conditions()
{
	if (Statesptr->startfile == ON)
	{
		read_depth_slopes();
		if (Statesptr->startq2d == ON) read_discharge_slopes();
	}
}

void dg2new::DG2Solver::read_depth_slopes()
{
	char h1x[256];
	char h1y[256];
	strcpy(h1x, Fnameptr->startfilename);
	strcat(h1x, "1x");
	strcpy(h1y, Fnameptr->startfilename);
	strcat(h1y, "1y");
	read_ascfile(h1x, Parptr, Arrptr->H1x, "Loading startfile 1x\n", verbose);
	read_ascfile(h1y, Parptr, Arrptr->H1y, "Loading startfile 1y\n", verbose);
}

void dg2new::DG2Solver::read_discharge_slopes()
{
	char hu1x[256];
	char hu1y[256];
	char hv1x[256];
	char hv1y[256];
	strcpy(hu1x, Fnameptr->startfilename);
	strcat(hu1x, ".Qx1x");
	strcpy(hu1y, Fnameptr->startfilename);
	strcat(hu1y, ".Qx1y");
	strcpy(hv1x, Fnameptr->startfilename);
	strcat(hv1x, ".Qy1x");
	strcpy(hv1y, Fnameptr->startfilename);
	strcat(hv1y, ".Qy1y");
	read_ascfile(hu1x, Parptr, Arrptr->HU1x, "Loading startfile Qx1x\n", verbose);
	read_ascfile(hu1y, Parptr, Arrptr->HU1y, "Loading startfile Qx1y\n", verbose);
	read_ascfile(hv1x, Parptr, Arrptr->HV1x, "Loading startfile Qy1x\n", verbose);
	read_ascfile(hv1y, Parptr, Arrptr->HV1y, "Loading startfile Qy1y\n", verbose);
}

void dg2new::DG2Solver::set_first_Tstep()
{
	if (Statesptr->adaptive_ts == ON)
	{
		Solverptr->Tstep = 1e-5;
	}
	else
	{
		Solverptr->Tstep = Solverptr->InitTstep;
	}
	Solverptr->MinTstep = Solverptr->InitTstep;
}

void dg2new::DG2Solver::print_Tstep()
{
    if (verbose == ON)
    {
        printf("t=%f\tdt=%f\n", Solverptr->t, Solverptr->Tstep);
    }
}

void dg2new::DG2Solver::apply_friction()
{
	if (Arrptr->Manningsn == nullptr && Parptr->FPn <= C(0.0)) return;

#pragma omp parallel for
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
            apply_friction(i, j);
        }
    }
}

void dg2new::DG2Solver::apply_friction(int i, int j)
{
	NUMERIC_TYPE H = Arrptr->H[j*Parptr->xsz+i];

    if (H <= Solverptr->DepthThresh) return;

    NUMERIC_TYPE n = (Arrptr->Manningsn == nullptr)
        ? Parptr->FPn : Arrptr->Manningsn[j*Parptr->xsz+i];

    NUMERIC_TYPE H1x = Arrptr->H1x[j*Parptr->xsz+i];
    NUMERIC_TYPE H1y = Arrptr->H1y[j*Parptr->xsz+i];
    NUMERIC_TYPE HU = Arrptr->HU[j*Parptr->xsz+i];
    NUMERIC_TYPE HV = Arrptr->HV[j*Parptr->xsz+i];
    NUMERIC_TYPE HU1x = Arrptr->HU1x[j*Parptr->xsz+i];
    NUMERIC_TYPE HU1y = Arrptr->HU1y[j*Parptr->xsz+i];
    NUMERIC_TYPE HV1x = Arrptr->HV1x[j*Parptr->xsz+i];
    NUMERIC_TYPE HV1y = Arrptr->HV1y[j*Parptr->xsz+i];

    NUMERIC_TYPE HU_friction = C(0.0);
    NUMERIC_TYPE HV_friction = C(0.0);
    NUMERIC_TYPE HU_Gx1_friction = C(0.0);
    NUMERIC_TYPE HV_Gx1_friction = C(0.0);
    NUMERIC_TYPE HU_Gx2_friction = C(0.0);
    NUMERIC_TYPE HV_Gx2_friction = C(0.0);
    NUMERIC_TYPE HU_Gy1_friction = C(0.0);
    NUMERIC_TYPE HV_Gy1_friction = C(0.0);
    NUMERIC_TYPE HU_Gy2_friction = C(0.0);
    NUMERIC_TYPE HV_Gy2_friction = C(0.0);

    apply_friction(n, H, HU, HV, HU_friction, HV_friction);
    apply_friction(n, H-H1x, HU-HU1x, HV-HV1x,
            HU_Gx1_friction, HV_Gx1_friction);
    apply_friction(n, H+H1x, HU+HU1x, HV+HV1x,
            HU_Gx2_friction, HV_Gx2_friction);
    apply_friction(n, H-H1y, HU-HU1y, HV-HV1y,
            HU_Gy1_friction, HV_Gy1_friction);
    apply_friction(n, H+H1y, HU+HU1y, HV+HV1y,
            HU_Gy2_friction, HV_Gy2_friction);

    Arrptr->HU[j*Parptr->xsz+i] = HU_friction;
    Arrptr->HV[j*Parptr->xsz+i] = HV_friction;
    Arrptr->HU1x[j*Parptr->xsz+i] = C(0.5)*(HU_Gx2_friction - HU_Gx1_friction);
    Arrptr->HV1x[j*Parptr->xsz+i] = C(0.5)*(HV_Gx2_friction - HV_Gx1_friction);
    Arrptr->HU1y[j*Parptr->xsz+i] = C(0.5)*(HU_Gy2_friction - HU_Gy1_friction);
    Arrptr->HV1y[j*Parptr->xsz+i] = C(0.5)*(HV_Gy2_friction - HV_Gy1_friction);
}

void dg2new::DG2Solver::apply_friction
(
    NUMERIC_TYPE n,
    NUMERIC_TYPE H,
    NUMERIC_TYPE HU_original,
    NUMERIC_TYPE HV_original,
    NUMERIC_TYPE& HU_frictional,
    NUMERIC_TYPE& HV_frictional
)
{
    if (H > Solverptr->DepthThresh)
    {
        NUMERIC_TYPE U = HU_original / H;
        NUMERIC_TYPE V = HV_original / H;

        if (FABS(U) <= Solverptr->SpeedThresh &&
                FABS(V) <= Solverptr->SpeedThresh)
        {
            HU_frictional = HU_original;
            HV_frictional = HV_original;
        }
        else
        {
            NUMERIC_TYPE Cf = Solverptr->g*n*n / pow(H, C(1.0)/C(3.0));

            NUMERIC_TYPE Sfx = - Cf*U*SQRT(U*U+V*V);
            NUMERIC_TYPE Sfy = - Cf*V*SQRT(U*U+V*V);

            NUMERIC_TYPE DDx = C(1.0) + Solverptr->Tstep*Cf /
                H*(C(2.0)*U*U + V*V)/SQRT(U*U+V*V);
            NUMERIC_TYPE DDy = C(1.0) + Solverptr->Tstep*Cf /
                H*(U*U + C(2.0)*V*V)/SQRT(U*U+V*V);

            HU_frictional = HU_original + Solverptr->Tstep*Sfx/DDx;
            HV_frictional = HV_original + Solverptr->Tstep*Sfy/DDy;
        }
    }
    else
    {
        HU_frictional = C(0.0);
        HV_frictional = C(0.0);
    }
}

// FIXME: BC treatment
void dg2new::DG2Solver::rk_stage1()
{
#pragma omp parallel for
	for (int j=1; j<Parptr->ysz-1; j++)
	{
		for(int i=1; i<Parptr->xsz-1; i++)
		{
            rk_stage1(i, j);
        }
    }
}

void dg2new::DG2Solver::rk_stage1(int i, int j)
{
    Increment L = space_operator(i, j, 0);

    Arrptr->H_int[j*Parptr->xsz+i] = Arrptr->H[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.H;
    Arrptr->HU_int[j*Parptr->xsz+i] = Arrptr->HU[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HU;
    Arrptr->HV_int[j*Parptr->xsz+i] = Arrptr->HV[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HV;

    Arrptr->H1x_int[j*Parptr->xsz+i] = Arrptr->H1x[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.H1x;
    Arrptr->HU1x_int[j*Parptr->xsz+i] = Arrptr->HU1x[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HU1x;
    Arrptr->HV1x_int[j*Parptr->xsz+i] = Arrptr->HV1x[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HV1x;

    Arrptr->H1y_int[j*Parptr->xsz+i] = Arrptr->H1y[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.H1y;
    Arrptr->HU1y_int[j*Parptr->xsz+i] = Arrptr->HU1y[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HU1y;
    Arrptr->HV1y_int[j*Parptr->xsz+i] = Arrptr->HV1y[j*Parptr->xsz+i]
        + Solverptr->Tstep*L.HV1y;

    if (Arrptr->H_int[j*Parptr->xsz+i] <= Solverptr->DepthThresh)
    {
        Arrptr->HU_int[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HU1x_int[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HU1y_int[j*Parptr->xsz+i] = C(0.0);

        Arrptr->HV_int[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HV1x_int[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HV1y_int[j*Parptr->xsz+i] = C(0.0);
    }
}

// FIXME: BC treatment
void dg2new::DG2Solver::rk_stage2()
{
#pragma omp parallel for
	for (int j=1; j<Parptr->ysz-1; j++)
	{
		for(int i=1; i<Parptr->xsz-1; i++)
		{
            rk_stage2(i, j);
        }
    }
}

void dg2new::DG2Solver::rk_stage2(int i, int j)
{
    Increment L = space_operator(i, j, 1);

    H_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->H[j*Parptr->xsz+i] +
            Arrptr->H_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.H);
    HU_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HU[j*Parptr->xsz+i] +
            Arrptr->HU_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HU);
    HV_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HV[j*Parptr->xsz+i] +
            Arrptr->HV_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HV);

    H1x_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->H1x[j*Parptr->xsz+i] +
            Arrptr->H1x_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.H1x);
    HU1x_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HU1x[j*Parptr->xsz+i] +
            Arrptr->HU1x_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HU1x);
    HV1x_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HV1x[j*Parptr->xsz+i] +
            Arrptr->HV1x_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HV1x);

    H1y_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->H1y[j*Parptr->xsz+i] +
            Arrptr->H1y_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.H1y);
    HU1y_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HU1y[j*Parptr->xsz+i] +
            Arrptr->HU1y_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HU1y);
    HV1y_new[j*Parptr->xsz+i] = C(0.5)*(Arrptr->HV1y[j*Parptr->xsz+i] +
            Arrptr->HV1y_int[j*Parptr->xsz+i] + Solverptr->Tstep*L.HV1y);

    if (H_new[j*Parptr->xsz+i] <= Solverptr->DepthThresh)
    {
        HU_new[j*Parptr->xsz+i] = C(0.0);
        HU1x_new[j*Parptr->xsz+i] = C(0.0);
        HU1y_new[j*Parptr->xsz+i] = C(0.0);

        HV_new[j*Parptr->xsz+i] = C(0.0);
        HV1x_new[j*Parptr->xsz+i] = C(0.0);
        HV1y_new[j*Parptr->xsz+i] = C(0.0);
    }
}

dg2new::Increment dg2new::DG2Solver::space_operator(int i, int j, int nmod)
{
    // FIXME: boundaries
    
    // east
    NUMERIC_TYPE DEM_east_LR = C(0.0);
    FlowVector star_east_L = {};
    FlowVector star_east_R = {};
    FlowVector F_east = {};
    {
        LocalFaceValue left = local_face_value(i, j, nmod, 2);
        LocalFaceValue right = local_face_value(i+1, j, nmod, 4);
        wetting_drying(left, right, 2, DEM_east_LR, star_east_L, star_east_R);

        HLL_x(Solverptr,
                star_east_L.H, star_east_L.HU, star_east_L.HV,
                star_east_R.H, star_east_R.HU, star_east_R.HV,
                F_east.H, F_east.HU, F_east.HV);
    }
    
    // west 
    NUMERIC_TYPE DEM_west_LR = C(0.0);
    FlowVector star_west_L = {};
    FlowVector star_west_R = {};
    FlowVector F_west = {};
    {
        LocalFaceValue right = local_face_value(i, j, nmod, 4);
        LocalFaceValue left = local_face_value(i-1, j, nmod, 2);
        wetting_drying(left, right, 4, DEM_west_LR, star_west_L, star_west_R);

        HLL_x(Solverptr,
                star_west_L.H, star_west_L.HU, star_west_L.HV,
                star_west_R.H, star_west_R.HU, star_west_R.HV,
                F_west.H, F_west.HU, F_west.HV);
    }

    // north
    NUMERIC_TYPE DEM_north_LR = C(0.0);
    FlowVector star_north_L = {};
    FlowVector star_north_R = {};
    FlowVector F_north = {};
    {
        LocalFaceValue left = local_face_value(i, j, nmod, 1);
        LocalFaceValue right = local_face_value(i, j-1, nmod, 3);
        wetting_drying(left, right, 1,
                DEM_north_LR, star_north_L, star_north_R);

        HLL_y(Solverptr,
                star_north_L.H, star_north_L.HU, star_north_L.HV,
                star_north_R.H, star_north_R.HU, star_north_R.HV,
                F_north.H, F_north.HU, F_north.HV);
    }
    
    // south
    NUMERIC_TYPE DEM_south_LR = C(0.0);
    FlowVector star_south_L = {};
    FlowVector star_south_R = {};
    FlowVector F_south = {};
    {
        LocalFaceValue right = local_face_value(i, j, nmod, 3);
        LocalFaceValue left = local_face_value(i, j+1, nmod, 1);
        wetting_drying(left, right, 3,
                DEM_south_LR, star_south_L, star_south_R);

        HLL_y(Solverptr,
                star_south_L.H, star_south_L.HU, star_south_L.HV,
                star_south_R.H, star_south_R.HU, star_south_R.HV,
                F_south.H, F_south.HU, F_south.HV);
    }
    
    // positivity preserving coefficients
    NUMERIC_TYPE DEM1x_bar = C(0.5)*(DEM_east_LR - DEM_west_LR)/SQRT(C(3.0));
    NUMERIC_TYPE DEM1y_bar = C(0.5)*(DEM_north_LR - DEM_south_LR)/SQRT(C(3.0));

    NUMERIC_TYPE H0x_bar = (star_east_L.H + star_west_R.H)/C(2.0);
    NUMERIC_TYPE H0y_bar = (star_north_L.H + star_south_R.H)/C(2.0);
    NUMERIC_TYPE H1x_bar = C(0.5)*(star_east_L.H - star_west_R.H)
        /SQRT(C(3.0));
    NUMERIC_TYPE H1y_bar = C(0.5)*(star_north_L.H - star_south_R.H)
        /SQRT(C(3.0));

    NUMERIC_TYPE HU0x_bar = (star_east_L.HU + star_west_R.HU)/C(2.0);
    NUMERIC_TYPE HU0y_bar = (star_north_L.HU + star_south_R.HU)/C(2.0);
    NUMERIC_TYPE HU1x_bar = C(0.5)*(star_east_L.HU - star_west_R.HU)
        /SQRT(C(3.0));
    NUMERIC_TYPE HU1y_bar = C(0.5)*(star_north_L.HU - star_south_R.HU)
        /SQRT(C(3.0));

    NUMERIC_TYPE HV0x_bar = (star_east_L.HV + star_west_R.HV)/C(2.0);
    NUMERIC_TYPE HV0y_bar = (star_north_L.HV + star_south_R.HV)/C(2.0);
    NUMERIC_TYPE HV1x_bar = C(0.5)*(star_east_L.HV - star_west_R.HV)
        /SQRT(C(3.0));
    NUMERIC_TYPE HV1y_bar = C(0.5)*(star_north_L.HV - star_south_R.HV)
        /SQRT(C(3.0));

    return dg2_operator(F_east, F_west, F_north, F_south,
            DEM1x_bar, DEM1y_bar,
            H0x_bar, H0y_bar, H1x_bar, H1y_bar,
            HU0x_bar, HU0y_bar, HU1x_bar, HU1y_bar,
            HV0x_bar, HV0y_bar, HV1x_bar, HV1y_bar);
}


dg2new::Increment dg2new::DG2Solver::dg2_operator
(
    FlowVector F_east,
    FlowVector F_west,
    FlowVector F_north,
    FlowVector F_south,
    NUMERIC_TYPE DEM1x_bar,
    NUMERIC_TYPE DEM1y_bar,
    NUMERIC_TYPE H0x_bar,
    NUMERIC_TYPE H0y_bar,
    NUMERIC_TYPE H1x_bar,
    NUMERIC_TYPE H1y_bar,
    NUMERIC_TYPE HU0x_bar,
    NUMERIC_TYPE HU0y_bar,
    NUMERIC_TYPE HU1x_bar,
    NUMERIC_TYPE HU1y_bar,
    NUMERIC_TYPE HV0x_bar,
    NUMERIC_TYPE HV0y_bar,
    NUMERIC_TYPE HV1x_bar,
    NUMERIC_TYPE HV1y_bar
)
{
    FlowVector SS0 =
    {
        C(0.0),
        - Solverptr->g*H0x_bar*C(2.0)*SQRT(C(3.0))*DEM1x_bar/Parptr->dx,
        - Solverptr->g*H0y_bar*C(2.0)*SQRT(C(3.0))*DEM1y_bar/Parptr->dy
    };
    FlowVector L0 = -(F_east - F_west)/Parptr->dx 
        - (F_north - F_south)/Parptr->dy + SS0;

    FlowVector F_Q1x =
        flux_x(H0x_bar-H1x_bar, HU0x_bar-HU1x_bar, HV0x_bar-HV1x_bar);
    FlowVector F_Q2x =
        flux_x(H0x_bar+H1x_bar, HU0x_bar+HU1x_bar, HV0x_bar+HV1x_bar);
    FlowVector SS1x =
    {
        C(0.0),
        C(2.0)*Solverptr->g*H1x_bar*DEM1x_bar,
        C(0.0)
    };
    FlowVector L1x = -SQRT(C(3.0))/Parptr->dx *
        (F_east + F_west - F_Q1x - F_Q2x + SS1x);

    FlowVector F_Q1y =
        flux_y(H0y_bar-H1y_bar, HU0y_bar-HU1y_bar, HV0y_bar-HV1y_bar);
    FlowVector F_Q2y =
        flux_y(H0y_bar+H1y_bar, HU0y_bar+HU1y_bar, HV0y_bar+HV1y_bar);
    FlowVector SS1y =
    {
        C(0.0),
        C(0.0),
        C(2.0)*Solverptr->g*H1y_bar*DEM1y_bar
    };
    FlowVector L1y = -SQRT(C(3.0))/Parptr->dy *
        (F_north + F_south - F_Q1y - F_Q2y + SS1y);

    return
    {
        L0.H, L1x.H, L1y.H,
        L0.HU, L1x.HU, L1y.HU,
        L0.HV, L1x.HV, L1y.HV
    };
}

dg2new::LocalFaceValue dg2new::DG2Solver::local_face_value
(
    int i,
    int j,
    int nmod,
    int ndir
)
{
    NUMERIC_TYPE DEM = Arrptr->DEM[j*Parptr->xsz+i];
    NUMERIC_TYPE DEM1x = Arrptr->DEM1x[j*Parptr->xsz+i];
    NUMERIC_TYPE DEM1y = Arrptr->DEM1y[j*Parptr->xsz+i];

    NUMERIC_TYPE H = C(0.0);
    NUMERIC_TYPE H1x = C(0.0);
    NUMERIC_TYPE H1y = C(0.0);
    NUMERIC_TYPE HU = C(0.0);
    NUMERIC_TYPE HU1x = C(0.0);
    NUMERIC_TYPE HU1y = C(0.0);
    NUMERIC_TYPE HV = C(0.0);
    NUMERIC_TYPE HV1x = C(0.0);
    NUMERIC_TYPE HV1y = C(0.0);
    
    if (nmod == 0)
    {
        H = Arrptr->H[j*Parptr->xsz+i];
        H1x = Arrptr->H1x[j*Parptr->xsz+i];
        H1y = Arrptr->H1y[j*Parptr->xsz+i];
        HU = Arrptr->HU[j*Parptr->xsz+i];
        HU1x = Arrptr->HU1x[j*Parptr->xsz+i];
        HU1y = Arrptr->HU1y[j*Parptr->xsz+i];
        HV = Arrptr->HV[j*Parptr->xsz+i];
        HV1x = Arrptr->HV1x[j*Parptr->xsz+i];
        HV1y = Arrptr->HV1y[j*Parptr->xsz+i];
    }
    else
    {
        H = Arrptr->H_int[j*Parptr->xsz+i];
        H1x = Arrptr->H1x_int[j*Parptr->xsz+i];
        H1y = Arrptr->H1y_int[j*Parptr->xsz+i];
        HU = Arrptr->HU_int[j*Parptr->xsz+i];
        HU1x = Arrptr->HU1x_int[j*Parptr->xsz+i];
        HU1y = Arrptr->HU1y_int[j*Parptr->xsz+i];
        HV = Arrptr->HV_int[j*Parptr->xsz+i];
        HV1x = Arrptr->HV1x_int[j*Parptr->xsz+i];
        HV1y = Arrptr->HV1y_int[j*Parptr->xsz+i];
    }

    NUMERIC_TYPE ETA = H + DEM;
    NUMERIC_TYPE ETA1x = H1x + DEM1x;
    NUMERIC_TYPE ETA1y = H1y + DEM1y;

    NUMERIC_TYPE H1_hat = C(0.0);
    NUMERIC_TYPE ETA1_hat = C(0.0);
    NUMERIC_TYPE HU1_hat = C(0.0);
    NUMERIC_TYPE HV1_hat = C(0.0);

    if (ndir == 2 || ndir == 4)
    {
        ETA1_hat = ETA1x;
        HU1_hat = HU1x;
        HV1_hat = HV1x;
        H1_hat = ETA1_hat - DEM1x;
    }
    else
    {
        ETA1_hat = ETA1y;
        HU1_hat = HU1y;
        HV1_hat = HV1y;
        H1_hat = ETA1_hat - DEM1y;
    }

    if (ndir == 1 || ndir == 2)
    {
        return
        {
            H + SQRT(C(3.0))*H1_hat,
            ETA + SQRT(C(3.0))*ETA1_hat,
            HU + SQRT(C(3.0))*HU1_hat,
            HV + SQRT(C(3.0))*HV1_hat
        };
    }
    else
    {
        return
        {
            H - SQRT(C(3.0))*H1_hat,
            ETA - SQRT(C(3.0))*ETA1_hat,
            HU - SQRT(C(3.0))*HU1_hat,
            HV - SQRT(C(3.0))*HV1_hat
        };
    }
}

void dg2new::DG2Solver::wetting_drying
(
    LocalFaceValue left,
    LocalFaceValue right,
    int ndir,
    NUMERIC_TYPE& DEM_LR,
    FlowVector& star_left,
    FlowVector& star_right
)
{
    NUMERIC_TYPE DEM_L = left.ETA - left.H;
    NUMERIC_TYPE DEM_R = right.ETA - right.H;
    
    NUMERIC_TYPE U_L = C(0.0);
    NUMERIC_TYPE V_L = C(0.0);
    if (left.H > Solverptr->DepthThresh)
    {
        U_L = left.HU / left.H;
        V_L = left.HV / left.H;
    }

    NUMERIC_TYPE U_R = C(0.0);
    NUMERIC_TYPE V_R = C(0.0);
    if (right.H > Solverptr->DepthThresh)
    {
        U_R = right.HU / right.H;
        V_R = right.HV / right.H;
    }

    DEM_LR = FMAX(DEM_L, DEM_R);

    NUMERIC_TYPE delta = C(0.0);
    if (ndir == 1 || ndir == 2)
    {
        delta = FMAX(C(0.0), -(left.ETA - DEM_LR));
    }
    else
    {
        delta = FMAX(C(0.0), -(right.ETA - DEM_LR));
    }

    star_left.H = FMAX(C(0.0), left.ETA - DEM_LR);
    star_left.HU = U_L * star_left.H;
    star_left.HV = V_L * star_left.H;

    star_right.H = FMAX(C(0.0), right.ETA - DEM_LR);
    star_right.HU = U_R * star_right.H;
    star_right.HV = V_R * star_right.H;
    
    if (delta > C(0.0))
    {
        DEM_LR -= delta;
    }
}

dg2new::FlowVector dg2new::DG2Solver::flux_x
(
    NUMERIC_TYPE H,
    NUMERIC_TYPE HU,
    NUMERIC_TYPE HV
)
{
    if (H > Solverptr->DepthThresh)
    {
        return
        {
            HU,
            HU*HU/H + (Solverptr->g/C(2.0))*(H*H),
            HU*HV/H
        };
    }
    else
    {
        return { C(0.0), C(0.0), C(0.0) };
    }
}

dg2new::FlowVector dg2new::DG2Solver::flux_y
(
    NUMERIC_TYPE H,
    NUMERIC_TYPE HU,
    NUMERIC_TYPE HV
)
{
    if (H > Solverptr->DepthThresh)
    {
        return
        {
            HU,
            HU*HV/H,
            HV*HV/H + (Solverptr->g/C(2.0))*(H*H)
        };
    }
    else
    {
        return { C(0.0), C(0.0), C(0.0) };
    }
}

void dg2new::DG2Solver::update_Tstep()
{
    if (Solverptr->t > C(0.0))
    {
        Solverptr->MinTstep = getmin(Solverptr->MinTstep, Solverptr->Tstep);
    }
    Solverptr->t += Solverptr->Tstep;
    Solverptr->itCount += 1;

    if (Statesptr->adaptive_ts == OFF) return; 

	NUMERIC_TYPE dt = Solverptr->InitTstep; // maximum timestep permitted

#ifndef _MSC_VER
#pragma omp parallel for reduction(min:dt)
#endif
	for (int j=0; j<Parptr->ysz; j++)
	{
		for(int i=0; i<Parptr->xsz; i++)
		{
            dt = update_Tstep(i, j, dt);
        }
    }

    Solverptr->Tstep = dt;
}

NUMERIC_TYPE dg2new::DG2Solver::update_Tstep(int i, int j, NUMERIC_TYPE dt)
{
    Arrptr->H[j*Parptr->xsz+i] = H_new[j*Parptr->xsz+i];
    Arrptr->HU[j*Parptr->xsz+i] = HU_new[j*Parptr->xsz+i];
    Arrptr->HV[j*Parptr->xsz+i] = HV_new[j*Parptr->xsz+i];

    Arrptr->H1x[j*Parptr->xsz+i] = H1x_new[j*Parptr->xsz+i];
    Arrptr->HU1x[j*Parptr->xsz+i] = HU1x_new[j*Parptr->xsz+i];
    Arrptr->HV1x[j*Parptr->xsz+i] = HV1x_new[j*Parptr->xsz+i];

    Arrptr->H1y[j*Parptr->xsz+i] = H1y_new[j*Parptr->xsz+i];
    Arrptr->HU1y[j*Parptr->xsz+i] = HU1y_new[j*Parptr->xsz+i];
    Arrptr->HV1y[j*Parptr->xsz+i] = HV1y_new[j*Parptr->xsz+i];

    if (Arrptr->H[j*Parptr->xsz+i] <= Solverptr->DepthThresh)
    {
        Arrptr->HU[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HU1x[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HU1y[j*Parptr->xsz+i] = C(0.0);

        Arrptr->HV[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HV1x[j*Parptr->xsz+i] = C(0.0);
        Arrptr->HV1y[j*Parptr->xsz+i] = C(0.0);

        return dt; // dry cells don't have a timestep
    }
    else
    {
        NUMERIC_TYPE U = Arrptr->HU[j*Parptr->xsz+i]/Arrptr->H[j*Parptr->xsz+i];
        NUMERIC_TYPE V = Arrptr->HV[j*Parptr->xsz+i]/Arrptr->H[j*Parptr->xsz+i];

        NUMERIC_TYPE dt_x = Solverptr->cfl*Parptr->dx /
            (FABS(U)+SQRT(Solverptr->g*Arrptr->H[j*Parptr->xsz+i]));
        NUMERIC_TYPE dt_y = Solverptr->cfl*Parptr->dy /
            (FABS(V)+SQRT(Solverptr->g*Arrptr->H[j*Parptr->xsz+i]));

        return std::min({dt, dt_x, dt_y});
    }

    return dt;
}

void dg2new::DG2Solver::update_mass_stats()
{
    if (Solverptr->t >= Parptr->MassTotal) {
        ::update_mass_stats(Statesptr, Parptr, Solverptr, BCptr, Arrptr);
        write_mass_stats(Fptr, Parptr, Solverptr, BCptr);
    }
}

void dg2new::DG2Solver::write_solution()
{
    if (Solverptr->t >= Parptr->SaveTotal) {
        if (Statesptr->voutput == ON)
        {
            update_velocity(Parptr, Solverptr, Arrptr);
        }
        write_solution_slopes(Fnameptr, Statesptr, Parptr, Solverptr, Arrptr);
        ::write_solution(Fnameptr, Statesptr, Parptr, Solverptr, Arrptr);
    }
}

void dg2new::DG2Solver::malloc_new_fields()
{
	H_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HU_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HV_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	H1x_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HU1x_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HV1x_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	H1y_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HU1y_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
	HV1y_new = memory_allocate_zero_numeric_legacy(Parptr->xsz*Parptr->ysz);
}
