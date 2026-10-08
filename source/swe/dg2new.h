#pragma once
#include "../lisflood.h"

namespace dg2new
{
    struct Increment
    {
        NUMERIC_TYPE H;
        NUMERIC_TYPE H1x;
        NUMERIC_TYPE H1y;
        NUMERIC_TYPE HU;
        NUMERIC_TYPE HU1x;
        NUMERIC_TYPE HU1y;
        NUMERIC_TYPE HV;
        NUMERIC_TYPE HV1x;
        NUMERIC_TYPE HV1y;
    };

    struct LocalFaceValue
    {
        NUMERIC_TYPE H;
        NUMERIC_TYPE ETA;
        NUMERIC_TYPE HU;
        NUMERIC_TYPE HV;
    };

    struct FlowVector
    {
        NUMERIC_TYPE H;
        NUMERIC_TYPE HU;
        NUMERIC_TYPE HV;
    };

    inline FlowVector operator-(const FlowVector& u)
    {
        return { -u.H, -u.HU, -u.HV };
    }
    
    inline FlowVector operator+(const FlowVector& lhs, const FlowVector& rhs)
    {
        return { lhs.H + rhs.H, lhs.HU + rhs.HU, lhs.HV + rhs.HV };
    }
    
    inline FlowVector operator-(const FlowVector& lhs, const FlowVector& rhs)
    {
        return { lhs.H - rhs.H, lhs.HU - rhs.HU, lhs.HV - rhs.HV };
    }
    
    inline FlowVector operator*(NUMERIC_TYPE scalar, const FlowVector& u)
    {
        return { scalar * u.H, scalar * u.HU, scalar * u.HV };
    }
    
    inline FlowVector operator/(const FlowVector& u, NUMERIC_TYPE scalar)
    {
        return { u.H / scalar, u.HU / scalar, u.HV / scalar };
    }

    class DG2Solver
    {
    public:
        DG2Solver
        (
            Fnames*,
            Files*,
            States*,
            Pars*,
            Solver*,
            BoundCs*,
            Stage*,
            Arrays*,
            int verbose
        );

        void solve();

        ~DG2Solver();

    private:
        void read_dem_slopes();
        void read_initial_conditions();
        void read_depth_slopes();
        void read_discharge_slopes();

        void set_first_Tstep();
        void print_Tstep();

        void apply_friction();
        void apply_friction(int i, int j);
        void apply_friction
        (
            NUMERIC_TYPE n,
            NUMERIC_TYPE H,
            NUMERIC_TYPE HU_original,
            NUMERIC_TYPE HV_original,
            NUMERIC_TYPE& HU_frictional,
            NUMERIC_TYPE& HV_frictional
        );

        void rk_stage1();
        void rk_stage1(int i, int j);
        void rk_stage2();
        void rk_stage2(int i, int j);

        /**
         * nmod == 0 denotes RK stage 1; 1 denotes RK stage 2
         */
        Increment space_operator(int i, int j, int nmod); 

        Increment dg2_operator
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
        );

        /**
         * nmod == 0: RK stage 1; 1: RK stage 2
         * ndir == 1: north face; 2: east face; 3: south face; 4: west face
         */
        LocalFaceValue local_face_value(int i, int j, int nmod, int ndir);

        void wetting_drying
        (
            LocalFaceValue left,
            LocalFaceValue right,
            int ndir,
            NUMERIC_TYPE& Z_LR,
            FlowVector& star_left,
            FlowVector& star_right
        );

        FlowVector flux_x(NUMERIC_TYPE H, NUMERIC_TYPE HU, NUMERIC_TYPE HV);
        FlowVector flux_y(NUMERIC_TYPE H, NUMERIC_TYPE HU, NUMERIC_TYPE HV);

        void update_Tstep();
        NUMERIC_TYPE update_Tstep(int i, int j, NUMERIC_TYPE dt);

        void update_mass_stats();
        void write_solution();

        void malloc_new_fields();

        Fnames* Fnameptr;
        Files* Fptr;
        States* Statesptr;
        Pars* Parptr;
        Solver* Solverptr;
        BoundCs* BCptr;
        Stage* Stageptr;
        Arrays* Arrptr;
        int verbose;

        NUMERIC_TYPE* H_new;
        NUMERIC_TYPE* HU_new;
        NUMERIC_TYPE* HV_new;
        NUMERIC_TYPE* H1x_new;
        NUMERIC_TYPE* HU1x_new;
        NUMERIC_TYPE* HV1x_new;
        NUMERIC_TYPE* H1y_new;
        NUMERIC_TYPE* HU1y_new;
        NUMERIC_TYPE* HV1y_new;
    };
}
