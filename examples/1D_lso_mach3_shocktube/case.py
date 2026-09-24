#!/usr/bin/env python3
import json
import math

# A right-moving Mach-3 shock initially located at x = 0.5.  The left state is
# the exact post-shock state for a stationary ideal-gas upstream state on the
# right; the discontinuity then provides a compact high-speed LSO smoke case.
gamma = 1.4
mach_shock = 3.0
rho_1 = 1.0
pres_1 = 1.0
sound_speed_1 = math.sqrt(gamma * pres_1 / rho_1)
rho_2 = rho_1 * (gamma + 1.0) * mach_shock**2 / ((gamma - 1.0) * mach_shock**2 + 2.0)
pres_2 = pres_1 * (1.0 + 2.0 * gamma / (gamma + 1.0) * (mach_shock**2 - 1.0))
vel_2 = mach_shock * sound_speed_1 * (1.0 - rho_1 / rho_2)

nx = 400
dx = 1.0 / nx
dt = 0.20 * dx / (mach_shock * sound_speed_1)

print(
    json.dumps(
        {
            # Logistics and domain
            "run_time_info": "T",
            "x_domain%beg": 0.0,
            "x_domain%end": 1.0,
            "m": nx - 1,
            "n": 0,
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": 60,
            "t_step_save": 20,
            # Shock-safe finite-difference setup
            "num_patches": 2,
            "model_eqns": "5eq",
            "num_fluids": 1,
            "alt_soundspeed": "F",
            "mpp_lim": "F",
            "mixture_err": "F",
            "time_stepper": "rk3",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "weno_avg": "T",
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "T",
            "riemann_solver": "hllc",
            "wave_speeds": "direct",
            "avg_state": "arithmetic",
            "fd_order": 4,
            "bc_x%beg": -3,
            "bc_x%end": -3,
            # LSO filter and 11-QOI output.  The shock-tube path has no IBM
            # radius; filter_sigma therefore defines the physical filter width.
            "lso_filter": "T",
            "lso_filter_wrt": "T",
            "lso_pp_filter": "T",
            "lso_stat_wrt": "T",
            "lso_closure_wrt": "T",
            "filter_sigma": 4.0 * dx,
            "lso_filter_sigma_target": 8.0 * dx,
            "lso_R_gas": 1.0,
            "lso_mu": 0.0,
            "fluid_pp(1)%k_therm": 0.0,
            # Output
            "format": "silo",
            "precision": "double",
            "prim_vars_wrt": "T",
            "E_wrt": "T",
            "parallel_io": "T",
            # Left/post-shock state
            "patch_icpp(1)%geometry": 1,
            "patch_icpp(1)%x_centroid": 0.25,
            "patch_icpp(1)%length_x": 0.5,
            "patch_icpp(1)%vel(1)": vel_2,
            "patch_icpp(1)%pres": pres_2,
            "patch_icpp(1)%alpha_rho(1)": rho_2,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Right/upstream state
            "patch_icpp(2)%geometry": 1,
            "patch_icpp(2)%x_centroid": 0.75,
            "patch_icpp(2)%length_x": 0.5,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%pres": pres_1,
            "patch_icpp(2)%alpha_rho(1)": rho_1,
            "patch_icpp(2)%alpha(1)": 1.0,
            # Ideal-gas convention: fluid_pp%gamma = 1 / (gamma_physical - 1).
            "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
            "fluid_pp(1)%eos": "ideal_gas",
            "fluid_pp(1)%cv": 1.0 / (gamma - 1.0),
        }
    )
)
