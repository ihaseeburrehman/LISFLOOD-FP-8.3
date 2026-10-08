//
// LISFLOOD-FP Metal FV1 Solver Kernels
// Ported from CUDA implementation
// Same numerical algorithm - different GPU language
//

#include <metal_stdlib>
using namespace metal;

// Constants matching LISFLOOD-FP definitions
#define DEPTH_THRESH 1e-3f
#define MIN_SAFE_H 1e-6f  // Minimum depth for safe division (1 micron)

// Portable atomic add for floats (works on any Metal version that supports atomic_uint)
inline void atomic_add_float(device atomic_uint* addr, float val) {
    uint expected = atomic_load_explicit(addr, memory_order_relaxed);
    while (true) {
        float f_expected = as_type<float>(expected);
        float f_new = f_expected + val;
        uint desired = as_type<uint>(f_new);
        if (atomic_compare_exchange_weak_explicit(addr, &expected, desired, memory_order_relaxed, memory_order_relaxed)) {
            break;
        }
    }
}

inline void account_boundary(device atomic_uint* totals, float outward_volume) {
    if (outward_volume > 0.0f) atomic_add_float(&totals[0], outward_volume);
    else if (outward_volume < 0.0f) atomic_add_float(&totals[1], -outward_volume);
}

// Boundary types matching ESourceType in lisflood.h
enum ENUM_BC_TYPE {
    BC_NONE = 0,   // closed / reflective
    BC_FREE = 1,   // transmissive
    BC_HFIX = 2,   // fixed water surface elevation
    BC_HVAR = 3,   // time-varying water surface elevation
    BC_QFIX = 4,   // fixed unit discharge (m2/s)
    BC_QVAR = 5,   // time-varying unit discharge (m2/s)
    BC_FREE6 = 6   // transmissive (alternate code)
};

// A boundary state split into the component normal to the edge and the one
// tangential to it, so the same routine serves both sweeps (the CPU does the
// same thing by swapping HU/HV in the argument list for the y boundaries).
struct EdgeState {
    float H;
    float Qn;  // momentum normal to the boundary
    float Qt;  // momentum tangential to the boundary
};

// Mirrors fv1::set_boundary_values() in swe/fv1.cpp.
//
// `in` is the interior cell's reconstructed state at the edge and `out` the
// state outside it; both may be modified, because the prescribed-discharge
// types also constrain the inside (the CPU does exactly this). `qn_sign` is +1
// where a positive prescribed discharge drives flow along this axis' positive
// direction and -1 where it opposes it, so that a positive BC value always
// means inflow on every edge.
inline void set_boundary_values(int type, float value, float Z,
                                float depth_thresh, float qn_sign,
                                thread EdgeState& in, thread EdgeState& out)
{
    if (type == BC_HFIX || type == BC_HVAR) {
        out.H = max(0.0f, value - Z);
        out.Qn = in.Qn;
        out.Qt = in.Qt;
    } else if (type == BC_QFIX || type == BC_QVAR) {
        out.Qn = qn_sign * value;
        if (abs(out.Qn) > 1e-12f) {
            in.H = max(in.H, 1.1f * depth_thresh);
            out.H = in.H;
        } else {
            out.H = in.H;
        }
        in.Qn = out.Qn;
        out.Qt = in.Qt;
    } else if (type == BC_FREE || type == BC_FREE6) {
        out.H = in.H;
        out.Qn = in.Qn;
        out.Qt = in.Qt;
    } else { // BC_NONE: closed wall
        out.H = in.H;
        out.Qn = -in.Qn;
        out.Qt = in.Qt;
    }
}

// Geometry parameters (set from host)
struct Geometry {
    int xsz;
    int ysz;
    int pitch;
    float dx;
    float dy;
    float tly;
    float blx;
    float bly;
};

// Physical parameters
struct PhysicalParams {
    float g;           // gravity
    float manning;     // Manning's n
    float speed_thresh;
};

// Solver parameters
struct SolverParams {
    float dt;
    float depth_thresh;
};

// Flow vector structure
struct FlowVector {
    float H;   // water height
    float HU;  // x-momentum
    float HV;  // y-momentum
};

// HLL Riemann solver in x-direction
inline FlowVector hll_x(FlowVector U_neg, FlowVector U_pos, float g, float depth_thresh) {
    float H_neg = max(U_neg.H, 0.0f);
    float H_pos = max(U_pos.H, 0.0f);
    
    if (H_neg <= depth_thresh && H_pos <= depth_thresh) {
        return FlowVector{0.0f, 0.0f, 0.0f};
    }

    float u_neg = (H_neg > depth_thresh) ? U_neg.HU / max(H_neg, 1e-6f) : 0.0f;
    float v_neg = (H_neg > depth_thresh) ? U_neg.HV / max(H_neg, 1e-6f) : 0.0f;
    float u_pos = (H_pos > depth_thresh) ? U_pos.HU / max(H_pos, 1e-6f) : 0.0f;
    float v_pos = (H_pos > depth_thresh) ? U_pos.HV / max(H_pos, 1e-6f) : 0.0f;
    
    float a_neg = sqrt(g * H_neg);
    float a_pos = sqrt(g * H_pos);
    
    // Toro estimate for H_star, U_star (matching hll_include.h)
    float h_star = pow(0.5f*(a_neg + a_pos) + 0.25f*(u_neg - u_pos), 2.0f) / g;
    float u_star = 0.5f*(u_neg + u_pos) + a_neg - a_pos;
    float a_star = sqrt(g * h_star);

    float s_neg, s_pos;
    if (H_neg <= depth_thresh) {
        s_neg = u_pos - 2.0f*a_pos;
    } else {
        s_neg = min(u_neg - a_neg, u_star - a_star);
    }

    if (H_pos <= depth_thresh) {
        s_pos = u_neg + 2.0f*a_neg;
    } else {
        s_pos = max(u_pos + a_pos, u_star + a_star);
    }
    
    s_neg = min(s_neg, 0.0f);
    s_pos = max(s_pos, 0.0f);
    
    FlowVector F_neg = {
        U_neg.HU,
        u_neg*U_neg.HU + 0.5f * g * H_neg * H_neg,
        U_neg.HV * u_neg
    };
    
    FlowVector F_pos = {
        U_pos.HU,
        u_pos*U_pos.HU + 0.5f * g * H_pos * H_pos,
        U_pos.HV * u_pos
    };
    
    float denom = s_pos - s_neg;
    if (abs(denom) < 1e-10f) return FlowVector{0.0f, 0.0f, 0.0f};
    
    FlowVector F;
    if (s_neg >= 0.0f) {
        F = F_neg;
    } else if (s_pos <= 0.0f) {
        F = F_pos;
    } else {
        F.H = (s_pos * F_neg.H - s_neg * F_pos.H + s_pos * s_neg * (U_pos.H - U_neg.H)) / denom;
        F.HU = (s_pos * F_neg.HU - s_neg * F_pos.HU + s_pos * s_neg * (U_pos.HU - U_neg.HU)) / denom;
        
        // Transverse momentum flux (HV) using S_mid logic from hll_include.h
        float s_mid = (s_neg * H_pos * (u_pos - s_pos) - s_pos * H_neg * (u_neg - s_neg)) 
                    / (H_pos * (u_pos - s_pos) - H_neg * (u_neg - s_neg));
        
        if (s_mid >= 0.0f) {
            F.HV = F.H * v_neg;
        } else {
            F.HV = F.H * v_pos;
        }
    }
    
    return F;
}

// HLL Riemann solver in y-direction
inline FlowVector hll_y(FlowVector U_neg, FlowVector U_pos, float g, float depth_thresh) {
    float H_neg = max(U_neg.H, 0.0f);
    float H_pos = max(U_pos.H, 0.0f);
    
    if (H_neg <= depth_thresh && H_pos <= depth_thresh) {
        return FlowVector{0.0f, 0.0f, 0.0f};
    }

    float u_neg = (H_neg > depth_thresh) ? U_neg.HU / max(H_neg, 1e-6f) : 0.0f;
    float v_neg = (H_neg > depth_thresh) ? U_neg.HV / max(H_neg, 1e-6f) : 0.0f;
    float u_pos = (H_pos > depth_thresh) ? U_pos.HU / max(H_pos, 1e-6f) : 0.0f;
    float v_pos = (H_pos > depth_thresh) ? U_pos.HV / max(H_pos, 1e-6f) : 0.0f;
    
    float a_neg = sqrt(g * H_neg);
    float a_pos = sqrt(g * H_pos);
    
    // Toro estimate for H_star, V_star (matching hll_include.h with V as primary)
    float h_star = pow(0.5f*(a_neg + a_pos) + 0.25f*(v_neg - v_pos), 2.0f) / g;
    float v_star = 0.5f*(v_neg + v_pos) + a_neg - a_pos;
    float a_star = sqrt(g * h_star);

    float s_neg, s_pos;
    if (H_neg <= depth_thresh) {
        s_neg = v_pos - 2.0f*a_pos;
    } else {
        s_neg = min(v_neg - a_neg, v_star - a_star);
    }

    if (H_pos <= depth_thresh) {
        s_pos = v_neg + 2.0f*a_neg;
    } else {
        s_pos = max(v_pos + a_pos, v_star + a_star);
    }
    
    s_neg = min(s_neg, 0.0f);
    s_pos = max(s_pos, 0.0f);
    
    FlowVector F_neg = {
        U_neg.HV,
        U_neg.HU * v_neg,
        v_neg * U_neg.HV + 0.5f * g * H_neg * H_neg
    };
    
    FlowVector F_pos = {
        U_pos.HV,
        U_pos.HU * v_pos,
        v_pos * U_pos.HV + 0.5f * g * H_pos * H_pos
    };
    
    float denom = s_pos - s_neg;
    if (abs(denom) < 1e-10f) return FlowVector{0.0f, 0.0f, 0.0f};
    
    FlowVector F;
    if (s_neg >= 0.0f) {
        F = F_neg;
    } else if (s_pos <= 0.0f) {
        F = F_pos;
    } else {
        F.H = (s_pos * F_neg.H - s_neg * F_pos.H + s_pos * s_neg * (U_pos.H - U_neg.H)) / denom;
        F.HV = (s_pos * F_neg.HV - s_neg * F_pos.HV + s_pos * s_neg * (U_pos.HV - U_neg.HV)) / denom;
        
        // Transverse momentum flux (HU in y-direction) using S_mid logic
        float s_mid = (s_neg * H_pos * (v_pos - s_pos) - s_pos * H_neg * (v_neg - s_neg)) 
                    / (H_pos * (v_pos - s_pos) - H_neg * (v_neg - s_neg));
        
        if (s_mid >= 0.0f) {
            F.HU = F.H * u_neg;
        } else {
            F.HU = F.H * u_pos;
        }
    }
    
    return F;
}

// Bed source term for x-direction
inline float bed_source_x(float Zstar_w, float Zstar_e, float Hstar_w, float Hstar_e, float ETA, float g, float dx) {
    float Zdagger_w = Zstar_w - max(0.0f, -(ETA - Zstar_w));
    float Zdagger_e = Zstar_e - max(0.0f, -(ETA - Zstar_e));
    return -g * 0.5f * (Hstar_w + Hstar_e) * (Zdagger_e - Zdagger_w) / dx;
}

// Bed source term for y-direction
inline float bed_source_y(float Zstar_n, float Zstar_s, float Hstar_n, float Hstar_s, float ETA, float g, float dy) {
    float Zdagger_n = Zstar_n - max(0.0f, -(ETA - Zstar_n));
    float Zdagger_s = Zstar_s - max(0.0f, -(ETA - Zstar_s));
    return -g * 0.5f * (Hstar_n + Hstar_s) * (Zdagger_s - Zdagger_n) / dy;
}

// Apply friction (Manning's equation)
kernel void apply_friction(
    device float* H [[buffer(0)]],
    device float* HU [[buffer(1)]],
    device float* HV [[buffer(2)]],
    constant Geometry& geom [[buffer(3)]],
    constant PhysicalParams& params [[buffer(4)]],
    constant SolverParams& solver [[buffer(5)]],
    device const float* manning [[buffer(6)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y + 1;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float h = H[idx];
    
    if (h > solver.depth_thresh) {
        float hu = HU[idx];
        float hv = HV[idx];
        
        // Use safe depth for division to prevent numerical instability
        float h_safe = max(h, MIN_SAFE_H);
        
        float n = (manning != nullptr) ? manning[idx] : params.manning;
        float u = hu / h_safe;
        float v = hv / h_safe;
        float speed = sqrt(u * u + v * v);
        float Cf = params.g * n * n / pow(h_safe, 0.33333333f);

        // Speed threshold - zero out momentum for very slow flow (matches fv1.cpp)
        if (abs(u) <= params.speed_thresh && abs(v) <= params.speed_thresh) {
            HU[idx] = 0.0f;
            HV[idx] = 0.0f;
            return;
        }
        
        // CPU semi-implicit friction (fv1.cpp lines 117-128)
        // The CPU uses: HU_new = HU_old + dt * Sf / D
        // where Sf = -Cf * u * speed and D = 1 + dt*Cf/h * directional_term
        float speed_safe = speed + 1e-12f; // Very small epsilon
        float denom_x = 1.0f + solver.dt * Cf / h_safe * (2.0f * u * u + v * v) / speed_safe;
        float denom_y = 1.0f + solver.dt * Cf / h_safe * (u * u + 2.0f * v * v) / speed_safe;
        
        float Sf_x = -Cf * u * speed;
        float Sf_y = -Cf * v * speed;
        
        HU[idx] = hu + solver.dt * Sf_x / denom_x;
        HV[idx] = hv + solver.dt * Sf_y / denom_y;
        
    } else {
        HU[idx] = 0.0f;
        HV[idx] = 0.0f;
    }
}

// Compute Z* at cell interfaces in x-direction
kernel void compute_zstar_x(
    device const float* DEM [[buffer(0)]],
    device float* Zstar_x [[buffer(1)]],
    constant Geometry& geom [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x;
    int j = gid.y + 1;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float Z_neg = DEM[idx];
    float Z_pos = DEM[idx + 1];
    
    Zstar_x[idx] = max(Z_neg, Z_pos);
}

// Compute Z* at cell interfaces in y-direction
kernel void compute_zstar_y(
    device const float* DEM [[buffer(0)]],
    device float* Zstar_y [[buffer(1)]],
    constant Geometry& geom [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float Z_neg = DEM[idx];         // North cell (row j)
    float Z_pos = DEM[(j+1) * geom.pitch + i]; // South cell (row j+1)
    
    Zstar_y[idx] = max(Z_neg, Z_pos);
}

// Threadgroup tile for the flux-sharing sweeps. Must match the host dispatch.
#define TG_X 16
#define TG_Y 16

// Reconstruct a cell's state at one interface (the "star" state).
inline FlowVector star_state(FlowVector U, float Z, float Zstar, float depth_thresh) {
    float u = (U.H > depth_thresh) ? U.HU / max(U.H, 1e-6f) : 0.0f;
    float v = (U.H > depth_thresh) ? U.HV / max(U.H, 1e-6f) : 0.0f;
    FlowVector S;
    S.H = max(0.0f, U.H + Z - Zstar);
    S.HU = S.H * u;
    S.HV = S.H * v;
    return S;
}

// Update flow variables in x-direction.
//
// Each thread evaluates exactly ONE interface flux (the interface between cell
// i and cell i+1) and publishes it to threadgroup memory; cell i then reads its
// west flux back as the *same bits* its neighbour computed. This mirrors
// cuda/fv1/cuda_fv1_solver.cu, and it is what makes the scheme conservative:
// evaluating the same interface independently in both adjacent cells is only
// algebraically equal, not bitwise equal, so the water leaving one cell did not
// match the water entering the next and the domain leaked mass.
// Threadgroups overlap by one column in x so every interface is covered;
// threads with tid.x == 0 contribute a flux but update no cell.
kernel void update_flow_x(
    device const float* H_old [[buffer(0)]],
    device const float* HU_old [[buffer(1)]],
    device const float* HV_old [[buffer(2)]],
    device float* H_new [[buffer(3)]],
    device float* HU_new [[buffer(4)]],
    device float* HV_new [[buffer(5)]],
    device const float* DEM [[buffer(6)]],
    device const float* Zstar_x [[buffer(7)]],
    constant Geometry& geom [[buffer(8)]],
    constant PhysicalParams& params [[buffer(9)]],
    constant SolverParams& solver [[buffer(10)]],
    device const int* boundary_types [[buffer(11)]],
    device const float* boundary_values [[buffer(12)]],
    device atomic_uint* mass_totals [[buffer(13)]],
    uint2 tgid [[threadgroup_position_in_grid]],
    uint2 tid [[thread_position_in_threadgroup]])
{
    threadgroup FlowVector F[TG_Y][TG_X];
    threadgroup float Hstar_sh[TG_Y][TG_X];

    int i = (int)(tgid.x * (TG_X - 1) + tid.x);   // interface / cell index 0..xsz
    int j = (int)(tgid.y * TG_Y + tid.y) + 1;     // row 1..ysz

    bool active = (i <= geom.xsz && j <= geom.ysz);

    FlowVector U_c, Ustar_c_e;
    float Z_c = 0.0f, Zstar_e = 0.0f;
    FlowVector F_e = FlowVector{0.0f, 0.0f, 0.0f};

    if (active) {
        int idx = j * geom.pitch + i;
        Zstar_e = Zstar_x[idx];

        Z_c = DEM[idx];
        U_c = FlowVector{H_old[idx], HU_old[idx], HV_old[idx]};
        Ustar_c_e = star_state(U_c, Z_c, Zstar_e, solver.depth_thresh);

        float Z_e = DEM[idx + 1];
        FlowVector U_e = FlowVector{H_old[idx + 1], HU_old[idx + 1], HV_old[idx + 1]};
        FlowVector Ustar_e = star_state(U_e, Z_e, Zstar_e, solver.depth_thresh);

        // Domain edges: replace the ghost-cell state with the prescribed
        // boundary condition. At i == 0 the outside state is Ustar_c_e (the
        // west ghost) and the interior is Ustar_e; at i == xsz it is the other
        // way round.
        if (i == 0 || i == geom.xsz) {
            bool west = (i == 0);
            int bc_idx = west ? (2 * geom.xsz + 2 * geom.ysz - (j - 1) - 1)
                              : (geom.xsz + (j - 1));
            int type = boundary_types[bc_idx];
            float value = boundary_values[bc_idx];
            float Z_in = west ? Z_e : Z_c;   // interior cell's bed elevation

            EdgeState in_s, out_s;
            if (west) {
                in_s = EdgeState{Ustar_e.H, Ustar_e.HU, Ustar_e.HV};
                out_s = EdgeState{Ustar_c_e.H, Ustar_c_e.HU, Ustar_c_e.HV};
            } else {
                in_s = EdgeState{Ustar_c_e.H, Ustar_c_e.HU, Ustar_c_e.HV};
                out_s = EdgeState{Ustar_e.H, Ustar_e.HU, Ustar_e.HV};
            }
            set_boundary_values(type, value, Z_in, solver.depth_thresh,
                                west ? 1.0f : -1.0f, in_s, out_s);
            if (west) {
                Ustar_e = FlowVector{in_s.H, in_s.Qn, in_s.Qt};
                Ustar_c_e = FlowVector{out_s.H, out_s.Qn, out_s.Qt};
            } else {
                Ustar_c_e = FlowVector{in_s.H, in_s.Qn, in_s.Qt};
                Ustar_e = FlowVector{out_s.H, out_s.Qn, out_s.Qt};
            }
        }

        F_e = hll_x(Ustar_c_e, Ustar_e, params.g, solver.depth_thresh);

        F[tid.y][tid.x] = F_e;
        Hstar_sh[tid.y][tid.x] = Ustar_e.H;

        // Domain-edge bookkeeping: i == 0 is the west face, i == xsz the east.
        if (i == 0) account_boundary(mass_totals, -F_e.H * solver.dt * geom.dy);
        if (i == geom.xsz) account_boundary(mass_totals, F_e.H * solver.dt * geom.dy);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid.x == 0 || !active || i < 1) return;

    int idx = j * geom.pitch + i;
    FlowVector F_w = F[tid.y][tid.x - 1];
    float Hstar_w = Hstar_sh[tid.y][tid.x - 1];
    float Zstar_w = Zstar_x[idx - 1];

    float ETA = U_c.H + Z_c;
    float bed_src = bed_source_x(Zstar_w, Zstar_e, Hstar_w, Ustar_c_e.H, ETA,
                                 params.g, geom.dx);

    // Partial (x-only) contribution of the unsplit FV1 update. No clamping and
    // no momentum zeroing here: update_flow_y still has to add the y-direction
    // divergence to this intermediate state.
    H_new[idx] = U_c.H - solver.dt * (F_e.H - F_w.H) / geom.dx;
    HU_new[idx] = U_c.HU - solver.dt * ((F_e.HU - F_w.HU) / geom.dx - bed_src);
    HV_new[idx] = U_c.HV - solver.dt * (F_e.HV - F_w.HV) / geom.dx;
}

// Update flow variables in y-direction (same flux-sharing scheme, overlapping
// in y). Metal's y axis points south (+j), so a thread at row j evaluates the
// interface between rows j and j+1 -- its own south face -- and cell j reads its
// north flux back from the thread above it.
kernel void update_flow_y(
    device const float* H_int [[buffer(0)]],
    device const float* HU_int [[buffer(1)]],
    device const float* HV_int [[buffer(2)]],
    device const float* H_old [[buffer(3)]],
    device const float* HU_old [[buffer(4)]],
    device const float* HV_old [[buffer(5)]],
    device float* H_new [[buffer(6)]],
    device float* HU_new [[buffer(7)]],
    device float* HV_new [[buffer(8)]],
    device const float* DEM [[buffer(9)]],
    device const float* Zstar_y [[buffer(10)]],
    constant Geometry& geom [[buffer(11)]],
    constant PhysicalParams& params [[buffer(12)]],
    constant SolverParams& solver [[buffer(13)]],
    device const int* boundary_types [[buffer(14)]],
    device const float* boundary_values [[buffer(15)]],
    device atomic_uint* mass_totals [[buffer(16)]],
    uint2 tgid [[threadgroup_position_in_grid]],
    uint2 tid [[thread_position_in_threadgroup]])
{
    threadgroup FlowVector F[TG_Y][TG_X];
    threadgroup float Hstar_sh[TG_Y][TG_X];

    int i = (int)(tgid.x * TG_X + tid.x) + 1;     // column 1..xsz
    int j = (int)(tgid.y * (TG_Y - 1) + tid.y);   // interface / cell row 0..ysz

    bool active = (i <= geom.xsz && j <= geom.ysz);

    FlowVector U_c, Ustar_c_s;
    float Z_c = 0.0f, Zstar_s = 0.0f;
    FlowVector F_s = FlowVector{0.0f, 0.0f, 0.0f};

    if (active) {
        int idx = j * geom.pitch + i;
        Zstar_s = Zstar_y[idx];

        Z_c = DEM[idx];
        U_c = FlowVector{H_old[idx], HU_old[idx], HV_old[idx]};
        Ustar_c_s = star_state(U_c, Z_c, Zstar_s, solver.depth_thresh);

        int idx_s = (j + 1) * geom.pitch + i;
        float Z_s = DEM[idx_s];
        FlowVector U_s = FlowVector{H_old[idx_s], HU_old[idx_s], HV_old[idx_s]};
        FlowVector Ustar_s = star_state(U_s, Z_s, Zstar_s, solver.depth_thresh);

        // Domain edges. Metal's y axis points south, so a positive prescribed
        // discharge is +Qn on the north edge and -Qn on the south edge; that
        // keeps "positive value == inflow" on every side, as on the CPU.
        if (j == 0 || j == geom.ysz) {
            bool north = (j == 0);
            int bc_idx = north ? (i - 1)
                               : (2 * geom.xsz + geom.ysz - (i - 1) - 1);
            int type = boundary_types[bc_idx];
            float value = boundary_values[bc_idx];
            float Z_in = north ? Z_s : Z_c;

            EdgeState in_s, out_s;
            if (north) {
                in_s = EdgeState{Ustar_s.H, Ustar_s.HV, Ustar_s.HU};
                out_s = EdgeState{Ustar_c_s.H, Ustar_c_s.HV, Ustar_c_s.HU};
            } else {
                in_s = EdgeState{Ustar_c_s.H, Ustar_c_s.HV, Ustar_c_s.HU};
                out_s = EdgeState{Ustar_s.H, Ustar_s.HV, Ustar_s.HU};
            }
            set_boundary_values(type, value, Z_in, solver.depth_thresh,
                                north ? 1.0f : -1.0f, in_s, out_s);
            if (north) {
                Ustar_s = FlowVector{in_s.H, in_s.Qt, in_s.Qn};
                Ustar_c_s = FlowVector{out_s.H, out_s.Qt, out_s.Qn};
            } else {
                Ustar_c_s = FlowVector{in_s.H, in_s.Qt, in_s.Qn};
                Ustar_s = FlowVector{out_s.H, out_s.Qt, out_s.Qn};
            }
        }

        F_s = hll_y(Ustar_c_s, Ustar_s, params.g, solver.depth_thresh);

        F[tid.y][tid.x] = F_s;
        Hstar_sh[tid.y][tid.x] = Ustar_s.H;

        if (j == 0) account_boundary(mass_totals, -F_s.H * solver.dt * geom.dx);
        if (j == geom.ysz) account_boundary(mass_totals, F_s.H * solver.dt * geom.dx);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid.y == 0 || !active || j < 1) return;

    int idx = j * geom.pitch + i;
    FlowVector F_n = F[tid.y - 1][tid.x];
    float Hstar_n = Hstar_sh[tid.y - 1][tid.x];
    float Zstar_n = Zstar_y[(j - 1) * geom.pitch + i];

    float ETA = U_c.H + Z_c;
    float bed_src = bed_source_y(Zstar_n, Zstar_s, Hstar_n, Ustar_c_s.H, ETA,
                                 params.g, geom.dy);

    float h_next = max(0.0f, H_int[idx] - solver.dt * (F_s.H - F_n.H) / geom.dy);
    float hu_next = HU_int[idx] - solver.dt * (F_s.HU - F_n.HU) / geom.dy;
    float hv_next = HV_int[idx] - solver.dt * ((F_s.HV - F_n.HV) / geom.dy - bed_src);

    // No "h <= depth_thresh => zero momentum" here: the CPU reference drops the
    // momentum of shallow cells in apply_friction, which runs after rainfall and
    // the point sources have topped the cell up.
    H_new[idx] = h_next;
    HU_new[idx] = hu_next;
    HV_new[idx] = hv_next;
}

// Compute CFL timestep per cell (for adaptive timestepping)
kernel void compute_dt_per_cell(
    device const float* H [[buffer(0)]],
    device const float* HU [[buffer(1)]],
    device const float* HV [[buffer(2)]],
    device float* dt_field [[buffer(3)]],
    constant Geometry& geom [[buffer(4)]],
    constant PhysicalParams& params [[buffer(5)]],
    constant float& cfl [[buffer(6)]],
    constant float& max_dt [[buffer(7)]],
    volatile device atomic_uint* min_dt_atomic [[buffer(8)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y + 1;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float h = H[idx];
    
    if (h > DEPTH_THRESH) {
        float hu = HU[idx];
        float hv = HV[idx];
        
        // Use safe depth for division to prevent numerical instability
        float h_safe = max(h, MIN_SAFE_H);
        float u = hu / h_safe;
        float v = hv / h_safe;
        float c = sqrt(params.g * h);
        
        float dt_x = geom.dx / (abs(u) + c);
        float dt_y = geom.dy / (abs(v) + c);
        
        dt_field[idx] = cfl * min(dt_x, dt_y);
    } else {
        dt_field[idx] = max_dt;
    }
    
    // Global reduction using SIMD and Atomics
    float local_dt = dt_field[idx];
    
    // SIMD reduction within warp/subgroup to reduce contention
    // Make sure to include metal_simdgroup if needed, but it's standard
    local_dt = simd_min(local_dt);
    
    if (simd_is_first()) {
        uint dt_uint = as_type<uint>(local_dt);
        atomic_fetch_min_explicit(min_dt_atomic, dt_uint, memory_order_relaxed);
    }
}

// Zero ghost cells (boundary)
kernel void zero_ghost_cells(
    device float* H [[buffer(0)]],
    device float* HU [[buffer(1)]],
    device float* HV [[buffer(2)]],
    constant Geometry& geom [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    // North and South boundaries
    if (gid < (uint)geom.pitch) {
        int i = gid;
        // North (j=0)
        H[i] = 0.0f;
        HU[i] = 0.0f;
        HV[i] = 0.0f;
        // South (j=ysz+1)
        int idx_s = (geom.ysz + 1) * geom.pitch + i;
        H[idx_s] = 0.0f;
        HU[idx_s] = 0.0f;
        HV[idx_s] = 0.0f;
    }
    
    // East and West boundaries
    if (gid < (uint)(geom.ysz + 2)) {
        int j = gid;
        // West (i=0)
        int idx_w = j * geom.pitch;
        H[idx_w] = 0.0f;
        HU[idx_w] = 0.0f;
        HV[idx_w] = 0.0f;
        // East (i=xsz+1)
        int idx_e = j * geom.pitch + geom.xsz + 1;
        H[idx_e] = 0.0f;
        HU[idx_e] = 0.0f;
        HV[idx_e] = 0.0f;
    }
}

// Apply uniform rainfall
kernel void apply_uniform_rain(
    device float* H [[buffer(0)]],
    device const float* DEM [[buffer(1)]],
    constant Geometry& geom [[buffer(2)]],
    constant float& rain_rate [[buffer(3)]],
    constant float& dt [[buffer(4)]],
    constant float& nodata [[buffer(5)]],
    device atomic_uint* mass_totals [[buffer(6)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y + 1;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float z = DEM[idx];
    
    // Don't add rain to nodata cells
    if (abs(z - nodata) > 1e-6f) {
        float increment = max(-H[idx], rain_rate * dt);
        H[idx] += increment;
        // Slot 3 is signed uniform-rain volume; negative rainfall is limited
        // by the water present, so accounting must use the applied increment.
        if (increment != 0.0f)
            atomic_add_float(&mass_totals[3], increment * geom.dx * geom.dy);
    }
}

// Helper for apply_rain (Dynamic)
float rate_at_cell(device const float* rain_data, 
                   int i, int j, 
                   constant Geometry& dem_geom, 
                   constant Geometry& rain_geom) 
{
    // i, j are 0-based DEM cell indices (inner domain)
    
    float x_offset = dem_geom.blx - rain_geom.blx;
    float tile_i_f = (x_offset + (float)i * dem_geom.dx) / rain_geom.dx;
    int tile_i = (int)floor(tile_i_f);
    
    float top_gap = dem_geom.tly - rain_geom.tly;
    float tile_j_f = ((float)j - top_gap/dem_geom.dy) * dem_geom.dy / rain_geom.dy;
    int tile_j = (int)floor(tile_j_f);
    
    if (tile_i >= 0 && tile_i < rain_geom.xsz && 
        tile_j >= 0 && tile_j < rain_geom.ysz) {
         // Rain data is row-major (y * xsz + x)
         return rain_data[tile_j * rain_geom.xsz + tile_i];
    }
    return 0.0f;
}

kernel void apply_rain(
    device float* H [[buffer(0)]],
    device const float* DEM [[buffer(1)]],
    device const float* rain_data [[buffer(2)]],
    constant Geometry& dem_geom [[buffer(3)]],
    constant Geometry& rain_geom [[buffer(4)]],
    constant float& dt [[buffer(5)]],
    constant float& nodata [[buffer(6)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x; // 0 to xsz-1
    int j = gid.y; // 0 to ysz-1
    
    if (i >= dem_geom.xsz || j >= dem_geom.ysz) return;
    
    // Calculate 1-based index for H/DEM which includes ghost cells
    int idx = (j + 1) * dem_geom.pitch + (i + 1);
    
    float z = DEM[idx];
    
    if (abs(z - nodata) < 1e-6f) return;
    
    float rate = rate_at_cell(rain_data, i, j, dem_geom, rain_geom);
    if (rate > 0.0f) {
        H[idx] += rate * dt;
    }
}

// Update ghost cells with reflective (wall) boundary conditions
// This copies interior values to ghost cells with velocity reflection
kernel void update_ghost_cells(
    device float* H [[buffer(0)]],
    device float* HU [[buffer(1)]],
    device float* HV [[buffer(2)]],
    constant Geometry& geom [[buffer(3)]],
    device const int* boundary_types [[buffer(4)]], // Packed array: N, E, S, W
    uint gid [[thread_position_in_grid]])
{
    int pitch = geom.pitch;
    int xsz = geom.xsz;
    int ysz = geom.ysz;
    
    // Offset for boundary_types packed array
    // North: 0
    // East: xsz
    // South: xsz + ysz
    // West: 2*xsz + ysz
    
    // West and East boundaries (i=0 and i=xsz+1)
    if (gid < (uint)(ysz + 2)) {
        int j = gid;
        
        // West boundary (i=0)
        if (j >= 1 && j <= ysz) {
            int bc_idx = 2 * xsz + ysz + (ysz - j); // LISFLOOD West index
            int type = boundary_types[bc_idx];
            
            int idx_ghost_w = j * pitch + 0;
            int idx_interior_w = j * pitch + 1;
            H[idx_ghost_w] = H[idx_interior_w];
            HV[idx_ghost_w] = HV[idx_interior_w];
            if (type == BC_FREE) {
                HU[idx_ghost_w] = HU[idx_interior_w]; // Simple transmissive (matches CPU fv1)
            } else {
                HU[idx_ghost_w] = -HU[idx_interior_w]; // Reflective
            }
        }
        
        // East boundary (i=xsz+1)
        if (j >= 1 && j <= ysz) {
            int bc_idx = xsz + (j - 1); // LISFLOOD East index
            int type = boundary_types[bc_idx];
            
            int idx_ghost_e = j * pitch + xsz + 1;
            int idx_interior_e = j * pitch + xsz;
            H[idx_ghost_e] = H[idx_interior_e];
            HV[idx_ghost_e] = HV[idx_interior_e];
            if (type == BC_FREE) {
                HU[idx_ghost_e] = HU[idx_interior_e]; // Simple transmissive (matches CPU fv1)
            } else {
                HU[idx_ghost_e] = -HU[idx_interior_e]; // Reflective
            }
        }
    }
    
    // North and South boundaries (j=0 and j=ysz+1)
    if (gid >= 1 && gid <= (uint)xsz) {
        int i = gid;
        
        // North boundary (j=0)
        int bc_idx_n = i - 1; // LISFLOOD North index
        int type_n = boundary_types[bc_idx_n];
        
        int idx_ghost_n = 0 * pitch + i;
        int idx_interior_n = 1 * pitch + i;
        H[idx_ghost_n] = H[idx_interior_n];
        HU[idx_ghost_n] = HU[idx_interior_n];
        if (type_n == BC_FREE) {
            HV[idx_ghost_n] = HV[idx_interior_n]; // Simple transmissive (matches CPU fv1)
        } else {
            HV[idx_ghost_n] = -HV[idx_interior_n]; // Reflective
        }
        
        // South boundary (j=ysz+1)
        int bc_idx_s = 2 * xsz + ysz - i; // LISFLOOD South index
        int type_s = boundary_types[bc_idx_s];
        
        int idx_ghost_s = (ysz + 1) * pitch + i;
        int idx_interior_s = ysz * pitch + i;
        H[idx_ghost_s] = H[idx_interior_s];
        HU[idx_ghost_s] = HU[idx_interior_s];
        if (type_s == BC_FREE) {
            HV[idx_ghost_s] = HV[idx_interior_s]; // Simple transmissive (matches CPU fv1)
        } else {
            HV[idx_ghost_s] = -HV[idx_interior_s]; // Reflective
        }
    }
}

// Drain water from NoData cells (implements drain_nodata functionality)
// Water entering NoData cells is eliminated (removed from simulation)
kernel void drain_nodata_cells(
    device float* H [[buffer(0)]],
    device float* HU [[buffer(1)]],
    device float* HV [[buffer(2)]],
    device const float* DEM [[buffer(3)]],
    constant Geometry& geom [[buffer(4)]],
    constant float& nodata [[buffer(5)]],
    device atomic_uint* mass_totals [[buffer(6)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y + 1;
    
    if (i > geom.xsz || j > geom.ysz) return;
    
    int idx = j * geom.pitch + i;
    float z = DEM[idx];
    
    // If this is a NoData cell, eliminate any water
    if (abs(z - nodata) < 1.0f) {
        if (H[idx] > 0.0f) atomic_add_float(&mass_totals[2], H[idx] * geom.dx * geom.dy);
        H[idx] = 0.0f;
        HU[idx] = 0.0f;
        HV[idx] = 0.0f;
    }
}

// Accumulate total outflow volume from FREE boundaries
// This kernel runs on 4 boundary sides and uses atomics to sum volume
kernel void accumulate_outflow(
    device const float* H [[buffer(0)]],
    device const float* HU [[buffer(1)]],
    device const float* HV [[buffer(2)]],
    device const int* boundary_types [[buffer(3)]],
    device atomic_uint* vol_out [[buffer(4)]], // Array of 4 totals: N, E, S, W
    constant Geometry& geom [[buffer(5)]],
    constant float& dt [[buffer(6)]],
    uint gid [[thread_position_in_grid]])
{
    int pitch = geom.pitch;
    int xsz = geom.xsz;
    int ysz = geom.ysz;
    
    // West and East boundaries
    if (gid < (uint)ysz) {
        int j = gid + 1;
        
        // West boundary (i=1) - outflow is NEGATIVE HU
        int bc_idx_w = 2 * xsz + ysz + (ysz - j);
        if (boundary_types[bc_idx_w] == BC_FREE) {
            float flux = -HU[j * pitch + 1] * geom.dy * dt;
            if (flux > 0) atomic_add_float(&vol_out[3], flux);
        }
        
        // East boundary (i=xsz) - outflow is POSITIVE HU
        int bc_idx_e = xsz + (j - 1);
        if (boundary_types[bc_idx_e] == BC_FREE) {
            float flux = HU[j * pitch + xsz] * geom.dy * dt;
            if (flux > 0) atomic_add_float(&vol_out[1], flux);
        }
    }
    
    // North and South boundaries
    if (gid < (uint)xsz) {
        int i = gid + 1;
        
        // North boundary (j=1) - outflow is NEGATIVE HV
        int bc_idx_n = i - 1;
        if (boundary_types[bc_idx_n] == BC_FREE) {
            float flux = -HV[1 * pitch + i] * geom.dx * dt;
            if (flux > 0) atomic_add_float(&vol_out[0], flux);
        }
        
        // South boundary (j=ysz) - outflow is POSITIVE HV
        int bc_idx_s = xsz + ysz + (xsz - i);
        if (boundary_types[bc_idx_s] == BC_FREE) {
            float flux = HV[ysz * pitch + i] * geom.dx * dt;
            if (flux > 0) atomic_add_float(&vol_out[2], flux);
        }
    }
}

// Track the running maximum water depth per cell.
// Mirrors update_max_field() in swe/stats.cpp, which the CPU solver calls once
// per timestep. Keeping this on the GPU avoids having to copy H back to the
// host every step just to maintain the .max output.
kernel void update_max_field(
    device const float* H [[buffer(0)]],
    device float* maxH [[buffer(1)]],
    constant Geometry& geom [[buffer(2)]],
    uint2 gid [[thread_position_in_grid]])
{
    int i = gid.x + 1;
    int j = gid.y + 1;

    if (i > geom.xsz || j > geom.ysz) return;

    int idx = j * geom.pitch + i;
    maxH[idx] = max(maxH[idx], H[idx]);
}
