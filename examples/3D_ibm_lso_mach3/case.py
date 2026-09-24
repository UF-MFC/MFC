import json
import math

# Stationary spherical particle in a uniform Mach-3 air stream.
# The domain contains 1000 stationary particles at 20% volume fraction. The
# resulting domain is intentionally large: at 50 cells/D it contains about
# 6.5e8 cells, so this is an HPC-scale case.
D = 0.1
num_particles = 1000
volume_fraction = 0.20
particle_volume = math.pi * D**3 / 6.0
domain_volume = num_particles * particle_volume / volume_fraction
lx = 25.0 * D
ly = math.sqrt(domain_volume / lx)
lz = ly
gamma = 1.4
rho_inf = 1.225
p_inf = 101325.0
c_inf = math.sqrt(gamma * p_inf / rho_inf)
mach_inf = 3.0
u_inf = mach_inf * c_inf

nx = round(lx / D * 50)
ny = round(ly / D * 50)
nz = round(lz / D * 50)
dx = D / 50.0
dt = 0.25 * dx / (u_inf + c_inf)


case = {
    # Logistics
    "run_time_info": "T",
    # Domain: 25D x Ly x Lz, with 50 cells/D.
    "x_domain%beg": -lx / 2.0,
    "x_domain%end": lx / 2.0,
    "y_domain%beg": -ly / 2.0,
    "y_domain%end": ly / 2.0,
    "z_domain%beg": -lz / 2.0,
    "z_domain%end": lz / 2.0,
    "cyl_coord": "F",
    "m": nx - 1,
    "n": ny - 1,
    "p": nz - 1,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": 5000,
    "t_step_save": 250,
    # Five-equation ideal-gas flow
    "num_patches": 1,
    "model_eqns": "5eq",
    "num_fluids": 1,
    "alt_soundspeed": "F",
    "mpp_lim": "F",
    "mixture_err": "T",
    "time_stepper": "rk3",
    # Shock-safe reconstruction
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": "arithmetic",
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": "hllc",
    "wave_speeds": "direct",
    # Uniform supersonic stream; the sphere generates the stationary bow shock.
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Stationary immersed sphere
    "ib": "T",
    "num_ibs": num_particles,
    "fd_order": 4,
    "viscous": "T",
    # LSO filtering and 11-QOI statistical products
    "lso_filter": "T",
    "lso_filter_wrt": "T",
    "lso_stat_wrt": "T",
    "filter_sigma": D / 4.0,
    "lso_R_gas": 287.05,
    "lso_mu": 1.84e-05,
    "fluid_pp(1)%k_therm": 0.0262,
    # Output
    "format": "silo",
    "precision": "double",
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "parallel_io": "T",
    # Uniform upstream air state
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": lx,
    "patch_icpp(1)%length_y": ly,
    "patch_icpp(1)%length_z": lz,
    "patch_icpp(1)%vel(1)": u_inf,
    "patch_icpp(1)%vel(2)": 0.0,
    "patch_icpp(1)%vel(3)": 0.0,
    "patch_icpp(1)%pres": p_inf,
    "patch_icpp(1)%alpha_rho(1)": rho_inf,
    "patch_icpp(1)%alpha(1)": 1.0,
    # Ideal-gas parameter convention: gamma parameter = 1/(gamma_physical - 1).
    "fluid_pp(1)%gamma": 1.0 / (gamma - 1.0),
    "fluid_pp(1)%eos": "ideal_gas",
    "fluid_pp(1)%cv": 287.05 / (gamma - 1.0),
    "fluid_pp(1)%Re(1)": rho_inf * u_inf * D / 1.84e-05,
}

# A deterministic 10 x 10 x 10 arrangement avoids overlaps while giving the
# requested particle count and volume fraction. The particles are stationary.
for ix in range(10):
    for iy in range(10):
        for iz in range(10):
            i = 1 + ix * 100 + iy * 10 + iz
            case[f"patch_ib({i})%geometry"] = 8
            case[f"patch_ib({i})%x_centroid"] = -lx / 2.0 + D / 2.0 + ix * (lx - D) / 9.0
            case[f"patch_ib({i})%y_centroid"] = -ly / 2.0 + D / 2.0 + iy * (ly - D) / 9.0
            case[f"patch_ib({i})%z_centroid"] = -lz / 2.0 + D / 2.0 + iz * (lz - D) / 9.0
            case[f"patch_ib({i})%radius"] = D / 2.0
            case[f"patch_ib({i})%vel(1)"] = 0.0
            case[f"patch_ib({i})%vel(2)"] = 0.0
            case[f"patch_ib({i})%vel(3)"] = 0.0
            case[f"patch_ib({i})%slip"] = "T"

print(json.dumps(case))
