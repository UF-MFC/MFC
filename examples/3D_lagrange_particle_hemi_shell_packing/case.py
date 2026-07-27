#!/usr/bin/env python3
import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument("--vf", type=float, default=0.10)
parser.add_argument("--seed", type=int, default=1)
parser.add_argument("--steps", type=int, default=1)
args, _ = parser.parse_known_args()
steps = max(1, args.steps)

R_INNER = 0.021
R_OUTER = 0.0225
DP = 100.0e-6

X_BEG = -0.025
X_END = 0.025
Y_BEG = -0.025
Y_END = 0.025
Z_BEG = 0.0
Z_END = 0.025

rho0 = 1.225
p0 = 101325.0
gamma = 1.4

data = {
    "run_time_info": "T",
    "x_domain%beg": X_BEG,
    "x_domain%end": X_END,
    "y_domain%beg": Y_BEG,
    "y_domain%end": Y_END,
    "z_domain%beg": Z_BEG,
    "z_domain%end": Z_END,
    "m": 63,
    "n": 63,
    "p": 31,
    "dt": 1.0e-8,
    "t_step_start": 0,
    "t_step_stop": steps,
    "t_step_save": steps,
    "model_eqns": 2,
    "alt_soundspeed": "F",
    "mixture_err": "F",
    "mpp_lim": "T",
    "time_stepper": 3,
    "weno_order": 5,
    "mapped_weno": "T",
    "mp_weno": "F",
    "avg_state": 2,
    "weno_eps": 1.0e-16,
    "riemann_solver": 2,
    "wave_speeds": 1,
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    "num_patches": 2,
    "num_fluids": 2,
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "parallel_io": "T",
    "lag_db_wrt": "T",
    "lag_id_wrt": "T",
    "lag_pos_wrt": "T",
    "lag_vel_wrt": "T",
    "lag_rad_wrt": "T",
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%pi_inf": 0.0,
    "fluid_pp(2)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(2)%pi_inf": 0.0,
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.5 * (X_BEG + X_END),
    "patch_icpp(1)%y_centroid": 0.5 * (Y_BEG + Y_END),
    "patch_icpp(1)%z_centroid": 0.5 * (Z_BEG + Z_END),
    "patch_icpp(1)%length_x": X_END - X_BEG,
    "patch_icpp(1)%length_y": Y_END - Y_BEG,
    "patch_icpp(1)%length_z": Z_END - Z_BEG,
    "patch_icpp(1)%vel(1)": 0.0,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p0,
    "patch_icpp(1)%alpha_rho(1)": rho0,
    "patch_icpp(1)%alpha_rho(2)": 1.0e-12,
    "patch_icpp(1)%alpha(1)": 1.0,
    "patch_icpp(1)%alpha(2)": 0.0,
    "patch_icpp(2)%geometry": 8,
    "patch_icpp(2)%alter_patch(1)": "T",
    "patch_icpp(2)%x_centroid": 0.0,
    "patch_icpp(2)%y_centroid": 0.0,
    "patch_icpp(2)%z_centroid": 0.0,
    "patch_icpp(2)%radius": R_INNER,
    "patch_icpp(2)%vel(1)": 0.0,
    "patch_icpp(2)%vel(2)": 0.0,
    "patch_icpp(2)%vel(3)": 0.0,
    "patch_icpp(2)%pres": 20.0 * p0,
    "patch_icpp(2)%alpha_rho(1)": rho0,
    "patch_icpp(2)%alpha_rho(2)": 1.0e-12,
    "patch_icpp(2)%alpha(1)": 1.0,
    "patch_icpp(2)%alpha(2)": 0.0,
    "particles_lagrange": "T",
    "fd_order": 4,
    "particle_pp%rho0ref_particle": 1740.0,
    "particle_pp%cp_particle": 1000.0,
    "particle_pp%ksp_col": 1.0,
    "particle_pp%nu_col": 0.3,
    "particle_pp%E_col": 1.0e9,
    "particle_pp%cor_col": 0.9,
    "lag_params%input_path": "input/lag_particles.dat",
    "lag_params%packing_flag": 2,
    "lag_params%packing_size_distribution": 0,
    "lag_params%packing_seed": args.seed,
    "lag_params%packing_volume_fraction": args.vf,
    "lag_params%packing_diameter_mean": DP,
    "lag_params%packing_min_spacing": 0.0,
    "lag_params%packing_centroid(1)": 0.0,
    "lag_params%packing_centroid(2)": 0.0,
    "lag_params%packing_centroid(3)": 0.0,
    "lag_params%packing_velocity(1)": 0.0,
    "lag_params%packing_velocity(2)": 0.0,
    "lag_params%packing_velocity(3)": 0.0,
    "lag_params%packing_shell_inner_radius": R_INNER,
    "lag_params%packing_shell_outer_radius": R_OUTER,
    "lag_params%vel_model": 2,
    "lag_params%drag_model": 1,
    "lag_params%solver_approach": 2,
    "lag_params%cluster_type": 1,
    "lag_params%pressure_corrector": "F",
    "lag_params%smooth_type": 1,
    "lag_params%pressure_force": "T",
    "lag_params%gravity_force": "F",
    "lag_params%qs_drag_model": 1,
    "lag_params%stokes_drag": 0,
    "lag_params%added_mass_model": 0,
    "lag_params%interpolation_order": 1,
    "lag_params%collision_force": "F",
    "lag_params%subcycle_collisions": "F",
    "lag_params%qs_fluct_force": "F",
    "lag_params%epsilonb": 1.0,
    "lag_params%valmaxvoid": 0.9,
    "lag_params%write_void_evol": "F",
    "lag_params%charwidth": Z_END - Z_BEG,
    "lag_params%charNz": 32,
}

print(json.dumps(data, indent=4))
