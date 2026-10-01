"""
MFC Parameter Definitions (Compact).

Single file containing all ~3,300 parameter definitions using loops.
"""

import re
from typing import Any, Dict

from .eos_families import EOS_FAMILIES
from .namelist_parser import get_fortran_constants
from .registry import REGISTRY, IndexedFamily
from .schema import ParamDef, ParamType

# Index limits — sourced from Fortran compile-time constants (m_constants.fpp).
# Falls back to the inline default when src/ is unavailable (e.g. Homebrew).
# Default must match src/common/m_constants.fpp — enforced by co-location.
_FC = get_fortran_constants()


def _fc(name: str, default: int) -> int:
    """Get a Fortran constant, using the inline default when m_constants.fpp is unavailable."""
    if _FC:
        if name not in _FC:
            raise RuntimeError(f"Fortran constant '{name}' not found in m_constants.fpp. Toolchain is out of sync with Fortran source.")
        return _FC[name]
    return default


NF = _fc("num_fluids_max", 10)  # fluid_pp
NPR = _fc("num_probes_max", 64)  # probe, acoustic
NB = _fc("num_bc_patches_max", 10)  # patch_bc
NUM_PATCHES_MAX = _fc("num_patches_max", 10)  # patch_icpp (Fortran array bound)
NIB = _fc("num_ib_patches_max_namelist", 54000)  # patch_ib namelist array bound
NAF = _fc("num_ib_airfoils_max", 5)  # ib_airfoil (Fortran array bound)
NSM = _fc("num_stl_models_max", 10)  # stl_models (Fortran array bound)
NPB = _fc("num_particle_clouds_max", 10)  # particle_cloud (Fortran array bound)
# Enumeration limits for families not yet converted to IndexedFamily.
# These are smaller than the Fortran array bounds to keep the registry compact.
# The CONSTRAINTS dict below uses the Fortran constants for validation.
NP = 10  # patch_icpp: has per-index variations, can't easily be IndexedFamily
NA = 4  # acoustic sources: enumerated individually


# Auto-generated Descriptions
# Descriptions are auto-generated from parameter names using naming conventions.
# Override with explicit desc= parameter when auto-generation is inadequate.

# Prefix descriptions for indexed parameter families
_PREFIX_DESCS = {
    "patch_icpp": "initial condition patch",
    "patch_ib": "immersed boundary",
    "patch_bc": "boundary condition patch",
    "fluid_pp": "fluid",
    "acoustic": "acoustic source",
    "probe": "probe",
    "integral": "integral region",
}

# Attribute descriptions (suffix after %)
_ATTR_DESCS = {
    # Geometry/position
    "geometry": "Geometry type",
    "x_centroid": "X-coordinate of centroid",
    "y_centroid": "Y-coordinate of centroid",
    "z_centroid": "Z-coordinate of centroid",
    "length_x": "X-dimension length",
    "length_y": "Y-dimension length",
    "length_z": "Z-dimension length",
    "radius": "Radius",
    "radii": "Radii array",
    "normal": "Normal direction",
    "theta": "Theta angle",
    "angles": "Orientation angles",
    # Physics
    "vel": "Velocity",
    "pres": "Pressure",
    "rho": "Density",
    "alpha": "Volume fraction",
    "alpha_rho": "Partial density",
    "gamma": "Specific heat ratio",
    "pi_inf": "Stiffness pressure",
    "cv": "Specific heat (const. volume)",
    "qv": "Heat of formation",
    "qvp": "Heat of formation prime",
    "G": "Shear modulus",
    "Re": "Reynolds number",
    "mul0": "Reference viscosity",
    "ss": "Surface tension",
    "pv": "Vapor pressure",
    # MHD
    "Bx": "Magnetic field (x-component)",
    "By": "Magnetic field (y-component)",
    "Bz": "Magnetic field (z-component)",
    # Model/smoothing
    "smoothen": "Enable smoothing",
    "smooth_patch_id": "Patch ID to smooth against",
    "smooth_coeff": "Smoothing coefficient",
    "alter_patch": "Alter with another patch",
    "model_filepath": "STL model file path",
    "model_spc": "Model spacing",
    "model_threshold": "Model threshold",
    "model_translate": "Model translation",
    "model_scale": "Model scale",
    "model_rotate": "Model rotation",
    # Bubbles
    "r0": "Initial bubble radius",
    "v0": "Initial bubble velocity",
    "p0": "Initial bubble pressure",
    "m0": "Initial bubble mass",
    # IB specific
    "slip": "Enable slip condition",
    "moving_ibm": "Enable moving boundary",
    "angular_vel": "Angular velocity",
    "mass": "Mass",
    # BC specific
    "vel_in": "Inlet velocity",
    "vel_out": "Outlet velocity",
    "alpha_rho_in": "Inlet partial density",
    "alpha_in": "Inlet volume fraction",
    "pres_in": "Inlet pressure",
    "pres_out": "Outlet pressure",
    "grcbc_in": "Enable GRCBC inlet",
    "grcbc_out": "Enable GRCBC outlet",
    "grcbc_vel_out": "Enable GRCBC velocity outlet",
    "isothermal_in": "Enable isothermal wall at the domain entrance (minimum coordinate)",
    "isothermal_out": "Enable isothermal wall at the domain exit (maximum coordinate)",
    # Acoustic
    "loc": "Location",
    "mag": "Magnitude",
    "pulse": "Pulse type",
    "support": "Support type",
    "frequency": "Frequency",
    "wavelength": "Wavelength",
    "length": "Length",
    "height": "Height",
    "delay": "Delay time",
    "dipole": "Enable dipole",
    "dir": "Direction",
    # Output
    "x": "X-coordinate",
    "y": "Y-coordinate",
    "z": "Z-coordinate",
    "xmin": "X minimum",
    "xmax": "X maximum",
    "ymin": "Y minimum",
    "ymax": "Y maximum",
    "zmin": "Z minimum",
    "zmax": "Z maximum",
    # Chemistry
    "Y": "Species mass fraction",
    # Shape coefficients
    "a": "Shape coefficient",
    # Elasticity
    "tau_e": "Elastic stress component",
    # Misc
    "cf_val": "Color function value",
    "hcid": "Hard-coded ID",
    "epsilon": "Interface thickness",
    "beta": "Shape parameter beta",
    "non_axis_sym": "Non-axisymmetric parameter",
}

# Simple parameter descriptions (non-indexed)
_SIMPLE_DESCS = {
    # Grid
    "m": "Grid cells in x-direction",
    "n": "Grid cells in y-direction",
    "p": "Grid cells in z-direction",
    "cyl_coord": "Enable cylindrical coordinates",
    "stretch_x": "Enable grid stretching in x",
    "stretch_y": "Enable grid stretching in y",
    "stretch_z": "Enable grid stretching in z",
    "a_x": "Grid stretching rate in x",
    "a_y": "Grid stretching rate in y",
    "a_z": "Grid stretching rate in z",
    "x_a": "Stretching start (negative x)",
    "x_b": "Stretching start (positive x)",
    "y_a": "Stretching start (negative y)",
    "y_b": "Stretching start (positive y)",
    "z_a": "Stretching start (negative z)",
    "z_b": "Stretching start (positive z)",
    "loops_x": "Stretching iterations in x",
    "loops_y": "Stretching iterations in y",
    "loops_z": "Stretching iterations in z",
    # Time
    "dt": "Time step size",
    "t_step_start": "Starting time step",
    "t_step_stop": "Ending time step",
    "t_step_save": "Save interval (steps)",
    "t_step_print": "Print interval (steps)",
    "t_stop": "Stop time",
    "t_save": "Save interval (time)",
    "time_stepper": "Time integration scheme",
    "cfl_target": "Target CFL number",
    "cfl_adap_dt": "Enable adaptive CFL time stepping",
    "cfl_const_dt": "Use constant CFL time stepping",
    "cfl_dt": "Enable CFL-based time stepping",
    "adap_dt": "Enable adaptive time stepping",
    "adap_dt_tol": "Adaptive time stepping tolerance",
    "adap_dt_max_iters": "Max iterations for adaptive dt",
    # Model
    "model_eqns": "Model equations",
    "num_fluids": "Number of fluids",
    "num_patches": "Number of IC patches",
    "mpp_lim": "Mixture pressure positivity limiter",
    # WENO
    "weno_order": "WENO reconstruction order",
    "weno_eps": "WENO epsilon parameter",
    "mapped_weno": "Enable mapped WENO",
    "wenoz": "Enable WENO-Z",
    "teno": "Enable TENO",
    "mp_weno": "Enable monotonicity-preserving WENO",
    # Riemann
    "riemann_solver": "Riemann solver",
    "wave_speeds": "Wave speed estimate method",
    "avg_state": "Average state",
    # Physics toggles
    "viscous": "Enable viscous effects",
    "mhd": "Enable magnetohydrodynamics",
    "hyper_cleaning": "Enable hyperbolic divergence cleaning",
    "hyper_cleaning_speed": "Divergence cleaning wave speed",
    "hyper_cleaning_tau": "Divergence cleaning damping time",
    "bubbles_euler": "Enable Euler bubble model",
    "bubbles_lagrange": "Enable Lagrangian bubbles",
    "polytropic": "Enable polytropic gas",
    "polydisperse": "Enable polydisperse bubbles",
    "qbmm": "Enable QBMM",
    "chemistry": "Enable chemistry",
    "surface_tension": "Enable surface tension",
    "hypoelasticity": "Enable hypoelastic model",
    "hyperelasticity": "Enable hyperelastic model",
    "relativity": "Enable special relativity",
    "ib": "Enable immersed boundaries",
    "collision_model": "Collision model for immersed boundaries (0=none, 1=soft sphere)",
    "coefficient_of_restitution": "Coefficient of restitution for IB collisions",
    "collision_time": "Characteristic collision time for IB collisions",
    "ib_coefficient_of_friction": "Coefficient of friction for IB collisions",
    # LSO variable-weight filter
    "lso_filter": "Enable LSO variable-weight Gaussian filter (applied at save steps)",
    "lso_filter_wrt": "Write LSO-filtered fields",
    "lso_down_sample_factor": "Coarsening stride for LSO-filtered output",
    "lso_stat_wrt": "Write LSO statistical products",
    "lso_R_gas": "Specific gas constant used by LSO statistical products",
    "lso_mu": "Dynamic viscosity used by LSO statistical products",
    "lso_filter_sigma_in": "Gaussian sigma already applied to post-process input",
    "lso_filter_sigma_target": "Target Gaussian sigma for post-process widening",
    "lso_pp2_n_passes_x": "Number of stage-2 post-process filter passes in x",
    "lso_pp2_n_passes_y": "Number of stage-2 post-process filter passes in y",
    "lso_pp2_n_passes_z": "Number of stage-2 post-process filter passes in z",
    "lso_pp2_a_x": "Stage-2 post-process per-pass stencil coefficients in x",
    "lso_pp2_a_y": "Stage-2 post-process per-pass stencil coefficients in y",
    "lso_pp2_a_z": "Stage-2 post-process per-pass stencil coefficients in z",
    "filter_sigma": "Target Gaussian filter standard deviation (physical units, same as domain coordinates)",
    "acoustic_source": "Enable acoustic sources",
    # Output
    "parallel_io": "Enable parallel I/O",
    "probe_wrt": "Write probe data",
    "prim_vars_wrt": "Write primitive variables",
    "cons_vars_wrt": "Write conservative variables",
    "run_time_info": "Print runtime info",
    "ib_state_wrt": "Write IB state and load data",
    # Misc
    "case_dir": "Case directory path",
    "cantera_file": "Cantera mechanism file",
    "num_ibs": "Number of immersed boundaries",
    "num_source": "Number of acoustic sources",
    "num_probes": "Number of probes",
    "num_integrals": "Number of integral regions",
    "nb": "Number of bubble bins",
    "R0ref": "Reference bubble radius",
    "sigma": "Surface tension coefficient",
    "Bx0": "Background magnetic field (x)",
    "old_grid": "Load grid from previous simulation",
    "old_ic": "Load initial conditions from previous",
    "t_step_old": "Time step to restart from",
    "fd_order": "Finite difference order",
    "recon_type": "Reconstruction type",
    "muscl_order": "MUSCL reconstruction order",
    "muscl_lim": "MUSCL limiter type",
    "muscl_eps": "MUSCL limiter slope-product threshold",
    "low_Mach": "Low Mach number correction",
    "bubble_model": "Bubble dynamics model",
    "Ca": "Cavitation number",
    "Web": "Weber number",
    "Re_inv": "Inverse Reynolds number",
    "format": "Output format",
    "precision": "Output precision",
    # Body forces
    "bf_x": "Enable body force in x",
    "bf_y": "Enable body force in y",
    "bf_z": "Enable body force in z",
    "k_x": "Body force wavenumber in x",
    "k_y": "Body force wavenumber in y",
    "k_z": "Body force wavenumber in z",
    "w_x": "Body force frequency in x",
    "w_y": "Body force frequency in y",
    "w_z": "Body force frequency in z",
    "p_x": "Body force phase in x",
    "p_y": "Body force phase in y",
    "p_z": "Body force phase in z",
    "g_x": "Gravitational acceleration in x",
    "g_y": "Gravitational acceleration in y",
    "g_z": "Gravitational acceleration in z",
    # More output
    "E_wrt": "Write energy field",
    "c_wrt": "Write sound speed field",
    "rho_wrt": "Write density field",
    "pres_wrt": "Write pressure field",
    "schlieren_wrt": "Write schlieren images",
    "cf_wrt": "Write color function",
    "omega_wrt": "Write vorticity",
    "qm_wrt": "Write Q-criterion",
    "liutex_wrt": "Write Liutex vortex field",
    "gamma_wrt": "Write gamma field",
    "heat_ratio_wrt": "Write heat capacity ratio",
    "pi_inf_wrt": "Write pi_inf field",
    "pres_inf_wrt": "Write reference pressure",
    "fft_wrt": "Write FFT output",
    "chem_wrt_T": "Write temperature (chemistry)",
    # Misc physics
    "alt_soundspeed": "Alternative sound speed formulation",
    "mixture_err": "Enable mixture error checking",
    "cont_damage": "Enable continuum damage model",
}


def _auto_describe(name: str) -> str:
    """Auto-generate description from parameter name."""
    # Check simple params first
    if name in _SIMPLE_DESCS:
        return _SIMPLE_DESCS[name]

    # Handle indexed params: prefix(N)%attr or prefix(N)%attr(M)
    match = re.match(r"([a-z_]+)\((\d+)\)%(.+)", name)
    if match:
        prefix, idx, attr = match.group(1), match.group(2), match.group(3)
        prefix_desc = _PREFIX_DESCS.get(prefix, prefix.replace("_", " "))

        # Check for nested index: attr(M) or attr(M, K)
        attr_match = re.match(r"([a-z_]+)\((\d+)(?:,\s*(\d+))?\)", attr)
        if attr_match:
            attr_base = attr_match.group(1)
            idx2 = attr_match.group(2)
            attr_desc = _ATTR_DESCS.get(attr_base, attr_base.replace("_", " "))
            return f"{attr_desc} {idx2} for {prefix_desc} {idx}"

        attr_desc = _ATTR_DESCS.get(attr, attr.replace("_", " "))
        return f"{attr_desc} for {prefix_desc} {idx}"

    # Handle bc_x%attr style (no index in prefix)
    if "%" in name:
        prefix, attr = name.split("%", 1)
        # Check for indexed attr
        attr_match = re.match(r"([a-z_]+)\((\d+)\)", attr)
        if attr_match:
            attr_base, idx = attr_match.group(1), attr_match.group(2)
            attr_desc = _ATTR_DESCS.get(attr_base, attr_base.replace("_", " "))
            return f"{attr_desc} {idx} for {prefix.replace('_', ' ')}"

        attr_desc = _ATTR_DESCS.get(attr, "")
        if attr_desc:
            return f"{attr_desc} for {prefix.replace('_', ' ')}"
        # Fallback: just clean up the name
        return f"{attr.replace('_', ' ').title()} for {prefix.replace('_', ' ')}"

    # Handle suffix-indexed: name(N) or name(N, M)
    match = re.match(r"([a-z_]+)\((\d+)(?:,\s*(\d+))?\)", name)
    if match:
        base, idx = match.group(1), match.group(2)
        # Handle _wrt patterns
        if base.endswith("_wrt"):
            field = base[:-4].replace("_", " ")
            return f"Write {field} for component {idx}"
        return f"{base.replace('_', ' ').title()} {idx}"

    # Fallback patterns
    if name.endswith("_wrt"):
        return f"Write {name[:-4].replace('_', ' ')}"
    if name.startswith("num_"):
        return f"Number of {name[4:].replace('_', ' ')}"

    # Last resort: clean up the name
    return name.replace("_", " ").replace("%", " ")


# Parameters that can be hard-coded for GPU case optimization
CASE_OPT_PARAMS = {
    "mapped_weno",
    "wenoz",
    "teno",
    "wenoz_q",
    "nb",
    "weno_order",
    "num_fluids",
    "mhd",
    "relativity",
    "igr_order",
    "viscous",
    "igr_iter_solver",
    "igr",
    "igr_pres_lim",
    "recon_type",
    "muscl_order",
    "muscl_lim",
}


# Data-driven Annotations for Doc Generation
# These dicts are the single source of truth for parameter hints in the docs.
# To annotate a new param, add an entry here instead of editing docs_gen.py.

HINTS = {
    "bc": {
        "vel_in_ramp": "Duration of the smooth start-up of the inflow velocity (0 = no ramp)",
        "vel_in_t0": "Time at which the inflow velocity ramp begins",
        "vel_in_frac0": "Fraction of the final inflow velocity held before the ramp",
        "grcbc_in": "Enables GRCBC subsonic inflow (bc type -7)",
        "grcbc_out": "Enables GRCBC subsonic outflow (bc type -8)",
        "grcbc_vel_out": "GRCBC velocity outlet (requires `grcbc_out`)",
        "vel_in": "Inlet velocity component (used with `grcbc_in`)",
        "vel_out": "Outlet velocity component (used with `grcbc_vel_out`)",
        "pres_in": "Inlet pressure (used with `grcbc_in`)",
        "pres_out": "Outlet pressure (used with `grcbc_out`)",
        "alpha_rho_in": "Inlet partial density per fluid (used with `grcbc_in`)",
        "alpha_in": "Inlet volume fraction per fluid (used with `grcbc_in`)",
        "vb1": "Boundary velocity component 1 at domain begin",
        "vb2": "Boundary velocity component 2 at domain begin",
        "vb3": "Boundary velocity component 3 at domain begin",
        "ve1": "Boundary velocity component 1 at domain end",
        "ve2": "Boundary velocity component 2 at domain end",
        "ve3": "Boundary velocity component 3 at domain end",
        "Twall_in": "Temperature of the entrance-side isothermal wall.",
        "Twall_out": "Temperature of the exit-side isothermal wall.",
    },
    "patch_bc": {
        "geometry": "Patch shape: 1=line, 2=circle, 3=rectangle",
        "type": "BC type applied within patch region",
        "dir": "Patch normal direction (1=x, 2=y, 3=z)",
        "loc": "Domain boundary (-1=begin, 1=end)",
        "centroid": "Patch center coordinate",
        "length": "Patch dimension",
        "radius": "Patch radius (geometry=2)",
    },
    "simplex_params": {
        "perturb_dens": "Enable simplex density perturbation",
        "perturb_dens_freq": "Density perturbation frequency",
        "perturb_dens_scale": "Density perturbation amplitude",
        "perturb_dens_offset": "Density perturbation offset seed",
        "perturb_vel": "Enable simplex velocity perturbation",
        "perturb_vel_freq": "Velocity perturbation frequency",
        "perturb_vel_scale": "Velocity perturbation amplitude",
        "perturb_vel_offset": "Velocity perturbation offset seed",
    },
    "fluid_pp": {
        "gamma": "Specific heat ratio (EOS)",
        "pi_inf": "Stiffness pressure (EOS)",
        "cv": "Specific heat at constant volume",
        "qv": "Heat of formation",
        "qvp": "Heat of formation derivative",
    },
}

# Tag → display name for docs. Dict order = priority when a param has multiple tags.
TAG_DISPLAY_NAMES = {
    "bubbles": "Bubble model",
    "mhd": "MHD",
    "chemistry": "Chemistry",
    "time": "Time-stepping",
    "grid": "Grid",
    "weno": "WENO",
    "viscosity": "Viscosity",
    "heat_conduction": "Heat conduction",
    "hypoelasticity": "Hypoelasticity",
    "surface_tension": "Surface tension",
    "acoustic": "Acoustic",
    "ib": "Immersed boundary",
    "reactive_burn": "Reactive burn",
    "probes": "Probe",
    "riemann": "Riemann solver",
    "relativity": "Relativity",
    "output": "Output",
    "bc": "Boundary condition",
}

# Prefix → hint for untagged simple params
PREFIX_HINTS = {
    "mixlayer_": "Mixing layer parameter",
    "nv_uvm_": "GPU memory management",
    "ic_": "Initial condition parameter",
}


def _lookup_hint(name):
    """Auto-derive constraint hint from HINTS dict using family+attribute matching."""
    if "%" not in name:
        # Check PREFIX_HINTS for simple params
        for prefix, label in PREFIX_HINTS.items():
            if name.startswith(prefix):
                return label
        return ""
    # Compound name: extract family and attribute
    prefix, attr_full = name.split("%", 1)
    # Normalize family: "bc_x" → "bc", "patch_bc(1)" → "patch"
    family = re.sub(r"[_(].*", "", prefix)
    if family not in HINTS:
        # Fallback: keep underscores — "patch_bc" → "patch_bc", "simplex_params" → "simplex_params"
        m = re.match(r"^[a-zA-Z_]+", prefix)
        family = m.group(0) if m else ""
    if family not in HINTS:
        return ""
    # Strip index from attr: "vel_in(1)" → "vel_in"
    m = re.match(r"^[a-zA-Z_0-9]+", attr_full)
    if not m:
        return ""
    attr = m.group(0)
    return HINTS[family].get(attr, "")


# Schema Validation for Constraints and Dependencies
# Uses rapidfuzz for "did you mean?" suggestions when typos are detected

_VALID_CONSTRAINT_KEYS = {"choices", "min", "max", "value_labels", "names", "fortran_prefix"}
_VALID_DEPENDENCY_KEYS = {"when_true", "when_set", "when_value"}
_VALID_CONDITION_KEYS = {"requires", "recommends", "requires_value"}


def _validate_constraint(param_name: str, constraint: Dict[str, Any]) -> None:
    """Validate a constraint dict has valid keys with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error

    invalid_keys = set(constraint.keys()) - _VALID_CONSTRAINT_KEYS
    if invalid_keys:
        # Get suggestion for the first invalid key
        first_invalid = next(iter(invalid_keys))
        raise ValueError(invalid_key_error(f"constraint for '{param_name}'", first_invalid, _VALID_CONSTRAINT_KEYS))

    # Validate types
    if "choices" in constraint and not isinstance(constraint["choices"], list):
        raise ValueError(f"Constraint 'choices' for '{param_name}' must be a list")
    if "min" in constraint and not isinstance(constraint["min"], (int, float)):
        raise ValueError(f"Constraint 'min' for '{param_name}' must be a number")
    if "max" in constraint and not isinstance(constraint["max"], (int, float)):
        raise ValueError(f"Constraint 'max' for '{param_name}' must be a number")
    if "value_labels" in constraint:
        if not isinstance(constraint["value_labels"], dict):
            raise ValueError(f"Constraint 'value_labels' for '{param_name}' must be a dict")
        if "choices" in constraint:
            for key in constraint["value_labels"]:
                if key not in constraint["choices"]:
                    raise ValueError(f"value_labels key {key!r} for '{param_name}' not in choices {constraint['choices']}")
    if "names" in constraint:
        names = constraint["names"]
        if not isinstance(names, dict):
            raise ValueError(f"Constraint 'names' for '{param_name}' must be a dict")
        for name, value in names.items():
            if not isinstance(name, str) or not re.match(r"^[a-z0-9][a-z0-9_]*$", name):
                raise ValueError(f"names key {name!r} for '{param_name}' must be a lowercase identifier")
            if not isinstance(value, int):
                raise ValueError(f"names value for '{param_name}'/{name!r} must be an int")
        if len(set(names.values())) != len(names):
            raise ValueError(f"names for '{param_name}' map two names to the same value")
        if "choices" in constraint and set(names.values()) != set(constraint["choices"]):
            raise ValueError(f"names for '{param_name}' must cover exactly its choices {constraint['choices']}")
    if "fortran_prefix" in constraint:
        if "names" not in constraint:
            raise ValueError(f"Constraint 'fortran_prefix' for '{param_name}' requires 'names'")
        if not isinstance(constraint["fortran_prefix"], str) or not re.match(r"^[a-z0-9][a-z0-9_]*$", constraint["fortran_prefix"]):
            raise ValueError(f"Constraint 'fortran_prefix' for '{param_name}' must be a lowercase identifier")


def _validate_dependency(param_name: str, dependency: Dict[str, Any]) -> None:
    """Validate a dependency dict has valid structure with 'did you mean?' suggestions."""
    # Import here to avoid circular import at module load time
    from .suggest import invalid_key_error

    invalid_keys = set(dependency.keys()) - _VALID_DEPENDENCY_KEYS
    if invalid_keys:
        first_invalid = next(iter(invalid_keys))
        raise ValueError(invalid_key_error(f"dependency for '{param_name}'", first_invalid, _VALID_DEPENDENCY_KEYS))

    def _validate_condition(cond_label: str, condition: Any) -> None:
        """Validate a condition dict (shared by when_true, when_set, when_value entries)."""
        if not isinstance(condition, dict):
            raise ValueError(f"Dependency '{cond_label}' for '{param_name}' must be a dict")
        invalid_cond_keys = set(condition.keys()) - _VALID_CONDITION_KEYS
        if invalid_cond_keys:
            first_invalid = next(iter(invalid_cond_keys))
            raise ValueError(invalid_key_error(f"condition in '{cond_label}' for '{param_name}'", first_invalid, _VALID_CONDITION_KEYS))
        for req_key in ["requires", "recommends"]:
            if req_key in condition and not isinstance(condition[req_key], list):
                raise ValueError(f"Dependency '{cond_label}/{req_key}' for '{param_name}' must be a list")
        if "requires_value" in condition:
            rv = condition["requires_value"]
            if not isinstance(rv, dict):
                raise ValueError(f"Dependency '{cond_label}/requires_value' for '{param_name}' must be a dict")
            for rv_param, rv_vals in rv.items():
                if not isinstance(rv_vals, list):
                    raise ValueError(f"Dependency '{cond_label}/requires_value/{rv_param}' for '{param_name}' must be a list")

    for condition_key in ["when_true", "when_set"]:
        if condition_key in dependency:
            _validate_condition(condition_key, dependency[condition_key])

    if "when_value" in dependency:
        wv = dependency["when_value"]
        if not isinstance(wv, dict):
            raise ValueError(f"Dependency 'when_value' for '{param_name}' must be a dict")
        for val, condition in wv.items():
            _validate_condition(f"when_value/{val}", condition)


def _validate_all_constraints(constraints: Dict[str, Dict]) -> None:
    """Validate all constraint definitions."""
    for param_name, constraint in constraints.items():
        _validate_constraint(param_name, constraint)


def _validate_all_dependencies(dependencies: Dict[str, Dict]) -> None:
    """Validate all dependency definitions."""
    for param_name, dependency in dependencies.items():
        _validate_dependency(param_name, dependency)


def get_value_label(param_name: str, value: int) -> str:
    """Look up the human-readable label for a parameter's integer code.

    Returns the label string, or ``str(value)`` when no label is defined.
    This is the single source of truth for value ↔ label mappings.
    """
    constraint = CONSTRAINTS.get(param_name)
    if constraint is None:
        return str(value)
    labels = constraint.get("value_labels")
    if labels is None:
        return str(value)
    return labels.get(value, str(value))


# Parameter constraints (choices, min, max)
CONSTRAINTS = {
    # Reconstruction
    "weno_order": {
        "choices": [0, 1, 3, 5, 7],
        "value_labels": {0: "MUSCL mode", 1: "1st order", 3: "WENO3", 5: "WENO5", 7: "WENO7"},
    },
    "recon_type": {
        "choices": [1, 2],
        "value_labels": {1: "WENO", 2: "MUSCL"},
        "names": {"weno": 1, "muscl": 2},
    },
    "muscl_order": {
        "choices": [1, 2],
        "value_labels": {1: "1st order", 2: "2nd order"},
        "names": {"first_order": 1, "second_order": 2},
    },
    "muscl_lim": {
        "choices": [0, 1, 2, 3, 4, 5],
        "value_labels": {0: "unlimited", 1: "minmod", 2: "MC", 3: "Van Albada", 4: "Van Leer", 5: "SUPERBEE"},
        "names": {"unlimited": 0, "minmod": 1, "mc": 2, "van_albada": 3, "van_leer": 4, "superbee": 5},
    },
    "int_comp": {
        "choices": [0, 1, 2],
        "value_labels": {0: "off", 1: "THINC", 2: "MTHINC"},
        "names": {"off": 0, "thinc": 1, "mthinc": 2},
    },
    # Time stepping
    "time_stepper": {
        "choices": [1, 2, 3],
        "value_labels": {1: "RK1 (Forward Euler)", 2: "RK2", 3: "RK3 (SSP)"},
        "names": {"rk1": 1, "rk2": 2, "rk3": 3},
    },
    # Riemann solver
    "riemann_solver": {
        "choices": [1, 2, 4, 5],
        "value_labels": {1: "HLL", 2: "HLLC", 4: "HLLD", 5: "Lax-Friedrichs"},
        "names": {"hll": 1, "hllc": 2, "hlld": 4, "lax_friedrichs": 5},
    },
    "wave_speeds": {
        "choices": [1, 2],
        "value_labels": {1: "direct", 2: "pressure"},
        "names": {"direct": 1, "pressure": 2},
    },
    "avg_state": {
        "choices": [1, 2],
        "value_labels": {1: "Roe", 2: "arithmetic"},
        "names": {"roe": 1, "arithmetic": 2},
    },
    # Model equations
    "model_eqns": {
        "choices": [1, 2, 3],
        "value_labels": {1: "Gamma-law", 2: "5-Equation", 3: "6-Equation"},
        "names": {"gamma_law": 1, "5eq": 2, "6eq": 3},
    },
    # Bubbles
    "bubble_model": {
        "choices": [0, 1, 2, 3],
        "value_labels": {0: "Particle", 1: "Gilmore", 2: "Keller-Miksis", 3: "Rayleigh-Plesset"},
        "names": {"particle": 0, "gilmore": 1, "keller_miksis": 2, "rayleigh_plesset": 3},
    },
    # Output
    "format": {
        "choices": [1, 2],
        "value_labels": {1: "Silo", 2: "binary"},
        "names": {"silo": 1, "binary": 2},
    },
    "precision": {
        "choices": [1, 2],
        "value_labels": {1: "single", 2: "double"},
        "names": {"single": 1, "double": 2},
    },
    # Time stepping (must be positive)
    "dt": {"min": 0},
    "t_stop": {"min": 0},
    "t_save": {"min": 0},
    "t_step_save": {"min": 1},
    "t_step_print": {"min": 1},
    "cfl_target": {"min": 0},
    "collision_temporal_resolution": {"min": 1},
    "ramp_ratio": {"min": 1},
    # WENO
    "weno_eps": {"min": 0},
    # MUSCL
    "muscl_eps": {"min": 0},
    # Physics (must be non-negative)
    "R0ref": {"min": 0},
    "sigma": {"min": 0},
    # Counts (must be positive)
    "num_fluids": {"min": 1, "max": NF},
    "num_patches": {"min": 0, "max": NUM_PATCHES_MAX},
    "num_ibs": {"min": 0},
    "ib_neighborhood_radius": {"min": 0},
    "num_source": {"min": 1},
    "num_probes": {"min": 1},
    "nb": {"min": 1},
    "m": {"min": 0},
    "n": {"min": 0},
    "p": {"min": 0},
}

# Parameter dependencies (requires, recommends)
DEPENDENCIES = {
    "bubbles_euler": {
        "when_true": {
            "recommends": ["nb", "polytropic"],
            "requires_value": {
                "model_eqns": [2],
                "riemann_solver": [2],
                "avg_state": [2],
            },
        }
    },
    "model_eqns": {
        "when_value": {
            2: {"requires": ["num_fluids"]},
            3: {"requires_value": {"riemann_solver": [2], "avg_state": [2], "wave_speeds": [1]}},
        }
    },
    "viscous": {
        "when_true": {
            "recommends": ["fluid_pp(1)%Re(1)"],
        }
    },
    "polydisperse": {
        "when_true": {
            "requires": ["nb", "poly_sigma"],
        }
    },
    "chemistry": {
        "when_true": {
            "requires": ["cantera_file"],
        }
    },
    "qbmm": {
        "when_true": {
            "recommends": ["bubbles_euler"],
        }
    },
    "ib": {
        "when_true": {
            "requires": ["num_ibs"],
        }
    },
    "collision_model": {
        "when_set": {
            "requires": ["ib", "coefficient_of_restitution", "collision_time"],
        }
    },
    "collision_temporal_resolution": {
        "when_set": {
            "requires": ["collision_model", "cfl_adap_dt"],
        }
    },
    "ramp_ratio": {
        "when_set": {
            "requires": ["cfl_adap_dt"],
        }
    },
    "acoustic_source": {
        "when_true": {
            "requires": ["num_source"],
        }
    },
    "probe_wrt": {
        "when_true": {
            "requires": ["num_probes", "fd_order"],
        }
    },
    "stretch_x": {
        "when_true": {
            "requires": ["a_x", "x_a", "x_b"],
        }
    },
    "stretch_y": {
        "when_true": {
            "requires": ["a_y", "y_a", "y_b"],
        }
    },
    "stretch_z": {
        "when_true": {
            "requires": ["a_z", "z_a", "z_b"],
        }
    },
    "bf_x": {
        "when_true": {
            "requires": ["k_x", "w_x", "p_x", "g_x"],
        }
    },
    "bf_y": {
        "when_true": {
            "requires": ["k_y", "w_y", "p_y", "g_y"],
        }
    },
    "bf_z": {
        "when_true": {
            "requires": ["k_z", "w_z", "p_z", "g_z"],
        }
    },
    "teno": {
        "when_true": {
            "requires": ["teno_CT"],
        }
    },
    "recon_type": {
        "when_value": {
            2: {"recommends": ["muscl_order", "muscl_lim"]},
        }
    },
    "surface_tension": {
        "when_true": {
            "requires": ["sigma"],
        }
    },
    "relativity": {
        "when_true": {
            "requires": ["mhd"],
        }
    },
    "riemann_hypo_ADC": {
        "when_true": {
            "requires": ["hypoelasticity"],
            "requires_value": {
                "riemann_solver": [2, 4],
            },
        }
    },
    "hypo_hll_interface_rhs": {
        "when_true": {
            "requires": ["hypoelasticity"],
            "requires_value": {
                "riemann_solver": [1],
            },
        }
    },
    "hll_u_interface": {
        "when_true": {
            "requires_value": {
                "riemann_solver": [1],
            },
        }
    },
    "schlieren_wrt": {
        "when_true": {
            "requires": ["fd_order"],
        }
    },
    "cfl_adap_dt": {
        "when_true": {
            "recommends": ["cfl_target"],
        }
    },
    "cfl_dt": {
        "when_true": {
            "recommends": ["cfl_target"],
        }
    },
}


def _r(name, ptype, tags=None, desc=None, hint=None, math=None, str_len=None, storage_precision=False):
    """Register a parameter with optional feature tags and description."""
    if hint is None:
        hint = _lookup_hint(name)
    if desc is None:
        from .descriptions import get_description

        desc = get_description(name)
    constraint = CONSTRAINTS.get(name)
    if constraint and "value_labels" in constraint:
        labels = constraint["value_labels"]
        by_value = {v: n for n, v in constraint.get("names", {}).items()}
        suffix = ", ".join(f"{v} ('{by_value[v]}')={labels[v]}" if v in by_value else f"{v}={labels[v]}" for v in sorted(labels))
        desc = f"{desc} ({suffix})".strip()
    REGISTRY.register(
        ParamDef(
            name=name,
            param_type=ptype,
            description=desc,
            case_optimization=(name in CASE_OPT_PARAMS),
            constraints=constraint,
            dependencies=DEPENDENCIES.get(name),
            tags=tags if tags else set(),
            hint=hint,
            math_symbol=math or "",
            str_len=str_len if str_len is not None else "name_len",
            storage_precision=storage_precision,
        )
    )


def _load():
    """Load all parameter definitions."""
    INT, REAL, LOG, STR = ParamType.INT, ParamType.REAL, ParamType.LOG, ParamType.STR
    A_REAL = ParamType.ANALYTIC_REAL

    # SIMPLE PARAMETERS (non-indexed)

    # Grid
    for n in ["m", "n", "p"]:
        _r(n, INT, {"grid"})
    _r("cyl_coord", LOG, {"grid"})
    for n in ["stretch_x", "stretch_y", "stretch_z"]:
        _r(n, LOG, {"grid"})
    for d in ["x", "y", "z"]:
        _r(f"{d}_a", REAL, {"grid"})
        _r(f"{d}_b", REAL, {"grid"})
        _r(f"a_{d}", REAL, {"grid"})
        _r(f"loops_{d}", INT, {"grid"})
        _r(f"{d}_domain%beg", REAL, {"grid"})
        _r(f"{d}_domain%end", REAL, {"grid"})

    # Time stepping
    for n in ["time_stepper", "t_step_old", "t_step_start", "t_step_stop", "t_step_save", "t_step_print", "adap_dt_max_iters"]:
        _r(n, INT, {"time"})
    _r("dt", REAL, {"time"}, math=r"\f$\Delta t\f$")
    _r("cfl_target", REAL, {"time"}, math=r"\f$\mathrm{CFL}\f$")
    for n in ["adap_dt_tol", "t_stop", "t_save", "ramp_ratio"]:
        _r(n, REAL, {"time"})
    for n in ["cfl_adap_dt", "cfl_const_dt", "cfl_dt", "adap_dt"]:
        _r(n, LOG, {"time"})

    # WENO/reconstruction
    _r("weno_order", INT, {"weno"})
    _r("recon_type", INT)
    _r("muscl_order", INT)
    _r("muscl_lim", INT)
    _r("muscl_eps", REAL)
    _r("weno_eps", REAL, {"weno"}, math=r"\f$\varepsilon\f$")
    _r("teno_CT", REAL, {"weno"}, math=r"\f$C_T\f$")
    _r("wenoz_q", REAL, {"weno"})
    for n in ["mapped_weno", "wenoz", "teno", "weno_avg", "mp_weno", "null_weights"]:
        _r(n, LOG, {"weno"})
    _r("weno_Re_flux", LOG, {"weno", "viscosity"})

    # Riemann solver
    for n in ["riemann_solver", "wave_speeds", "avg_state", "low_Mach"]:
        _r(n, INT, {"riemann"})

    # MHD
    _r("Bx0", REAL, {"mhd"}, math=r"\f$B_{x,0}\f$")
    _r("hyper_cleaning_speed", REAL, {"mhd"}, math=r"\f$c_h\f$")
    _r("hyper_cleaning_tau", REAL, {"mhd"})
    for n in ["mhd", "hyper_cleaning"]:
        _r(n, LOG, {"mhd"})

    # Bubbles
    _r("R0ref", REAL, {"bubbles"}, math=r"\f$R_0\f$")
    _r("nb", INT, {"bubbles"}, math=r"\f$N_b\f$")
    _r("Web", REAL, {"bubbles"}, math=r"\f$\mathrm{We}\f$")
    _r("Ca", REAL, {"bubbles"}, math=r"\f$\mathrm{Ca}\f$")
    _r("Re_inv", REAL, {"bubbles", "viscosity"}, math=r"\f$\mathrm{Re}^{-1}\f$")
    _r("bubble_model", INT, {"bubbles"})
    for n in ["polytropic", "bubbles_euler", "polydisperse", "qbmm", "bubbles_lagrange"]:
        _r(n, LOG, {"bubbles"})

    # Subgrid solid particles (Euler-Lagrange)
    _r("particles_lagrange", LOG, {"particles"})
    for a in ["rho0ref_particle", "cp_particle", "ksp_col", "nu_col", "E_col", "cor_col"]:
        _r(f"particle_pp%{a}", REAL, {"particles"})

    # Viscosity
    _r("viscous", LOG, {"viscosity"})

    # Hypoelasticity
    _r("hypoelasticity", LOG, {"hypoelasticity"})
    _r("riemann_hypo_ADC", LOG, {"hypoelasticity"})
    _r("ADC_kappa", REAL, {"hypoelasticity"})
    _r("hypo_hll_interface_rhs", LOG, {"hypoelasticity"})
    _r("hll_u_interface", LOG, {"riemann"})

    # Surface tension
    _r("sigma", REAL, {"surface_tension"}, math=r"\f$\sigma\f$")
    _r("surface_tension", LOG, {"surface_tension"})

    # JWL reaction and diagnostic controls
    _r("jwl_wrt", LOG, desc="Write JWL temperature, product fraction, and reaction progress")
    _r("jwl_afterburn", LOG, desc="Enable JWL afterburn energy release")
    _r("jwl_ab_model", INT, desc="JWL afterburn rate model")
    for n in ["jwl_q_ab", "jwl_ab_tau", "jwl_ab_A", "jwl_ab_theta", "jwl_ab_n"]:
        _r(n, REAL)
    _r("jwl_reactive", LOG, desc="Enable JWL++ pressure-driven reactive burn")
    for n in ["jwl_G", "jwl_b_exp"]:
        _r(n, REAL)
    _r("prog_burn", LOG, desc="Enable kinematic JWL program burn")
    for n in ["pb_D_cj", "pb_width", "pb_x_det", "pb_y_det", "pb_z_det", "pb_t_det"]:
        _r(n, REAL)

    # Chemistry
    _r("cantera_file", STR, {"chemistry"})
    _r("chemistry", LOG, {"chemistry"})

    # Condensed-phase reactive burn (programmed pressure burn on the multi-fluid model)
    _r("reactive_burn", LOG, {"reactive_burn"})
    _r("rburn%model", INT, {"reactive_burn"})
    for a in ["k", "pign", "pref", "n", "ta"]:
        _r(f"rburn%{a}", REAL, {"reactive_burn"})
    for a in ["rho0", "q", "ki", "kg", "m1", "m2", "n1", "n2", "n3"]:
        _r(f"rburn%{a}", REAL, {"reactive_burn"})
    _r("rburn%substeps", INT, {"reactive_burn"})

    # Acoustic
    _r("num_source", INT, {"acoustic"})
    _r("acoustic_source", LOG, {"acoustic"})

    # Immersed boundary
    _r("num_ibs", INT, {"ib"})
    _r("num_stl_models", INT, {"ib"})
    _r("num_particle_clouds", INT, {"ib"})
    _r("ib_neighborhood_radius", INT, {"ib"})
    _r("ib", LOG, {"ib"})
    _r("collision_model", INT, {"ib"})
    _r("collision_temporal_resolution", INT, {"ib"})
    _r("coefficient_of_restitution", REAL, {"ib"})
    _r("collision_time", REAL, {"ib"})
    _r("ib_coefficient_of_friction", REAL, {"ib"})
    _r("many_ib_patch_parallelism", LOG, {"ib"})

    # Probes
    _r("num_probes", INT, {"probes"})
    _r("probe_wrt", LOG, {"output", "probes"})

    # Output
    _r("precision", INT, {"output"})
    _r("format", INT, {"output"})
    _r("ib_force_stride", INT, {"output", "ib"})
    for n in ["parallel_io", "file_per_process", "run_time_info", "prim_vars_wrt", "cons_vars_wrt", "fft_wrt", "ib_state_wrt", "ib_force_wrt"]:
        _r(n, LOG, {"output"})

    # LSO variable-weight filter
    _r("lso_filter", LOG, {"filter"})
    _r("lso_filter_wrt", LOG, {"filter"})
    _r("filter_sigma", REAL, {"filter"})
    _r("lso_filter_sigma_in", REAL, {"filter"})
    _r("lso_filter_sigma_target", REAL, {"filter"})
    _r("lso_down_sample_factor", INT, {"filter"})
    _r("lso_stat_wrt", LOG, {"filter"})
    _r("lso_R_gas", REAL, {"filter"})
    _r("lso_mu", REAL, {"filter"})
    for n in ["lso_n_passes_x", "lso_n_passes_y", "lso_n_passes_z"]:
        _r(n, INT, {"filter"})
    for n in ["lso_a_x", "lso_a_y", "lso_a_z"]:
        _r(f"{n}(1)", REAL, {"filter"})
    for n in ["lso2_n_passes_x", "lso2_n_passes_y", "lso2_n_passes_z"]:
        _r(n, INT, {"filter"})
    for n in ["lso2_a_x", "lso2_a_y", "lso2_a_z"]:
        _r(f"{n}(1)", REAL, {"filter"})
    _r("lso_pp_filter", LOG, {"filter"})
    _r("lso_closure_wrt", LOG, {"filter"})
    for n in ["lso_pp_n_passes_x", "lso_pp_n_passes_y", "lso_pp_n_passes_z"]:
        _r(n, INT, {"filter"})
    for n in ["lso_pp_a_x", "lso_pp_a_y", "lso_pp_a_z"]:
        _r(f"{n}(1)", REAL, {"filter"})
    for n in ["lso_pp2_n_passes_x", "lso_pp2_n_passes_y", "lso_pp2_n_passes_z"]:
        _r(n, INT, {"filter"})
    for n in ["lso_pp2_a_x", "lso_pp2_a_y", "lso_pp2_a_z"]:
        _r(f"{n}(1)", REAL, {"filter"})
    for n in [
        "schlieren_wrt",
        "alpha_wrt",
        "rho_wrt",
        "E_wrt",
        "pres_wrt",
        "gamma_wrt",
        "heat_ratio_wrt",
        "pi_inf_wrt",
        "pres_inf_wrt",
        "c_wrt",
        "T_wrt",
        "qm_wrt",
        "liutex_wrt",
        "cf_wrt",
        "sim_data",
        "output_partial_domain",
    ]:
        _r(n, LOG, {"output"})
    for d in ["x", "y", "z"]:
        _r(f"{d}_output%beg", REAL, {"output"})
        _r(f"{d}_output%end", REAL, {"output"})
    # Lagrangian output
    for v in [
        "lag_header",
        "lag_txt_wrt",
        "lag_db_wrt",
        "lag_id_wrt",
        "lag_pos_wrt",
        "lag_pos_prev_wrt",
        "lag_vel_wrt",
        "lag_rad_wrt",
        "lag_rvel_wrt",
        "lag_r0_wrt",
        "lag_rmax_wrt",
        "lag_rmin_wrt",
        "lag_dphidt_wrt",
        "lag_pres_wrt",
        "lag_mv_wrt",
        "lag_mg_wrt",
        "lag_betaT_wrt",
        "lag_betaC_wrt",
    ]:
        _r(v, LOG, {"bubbles", "output"})

    # Boundary conditions
    for d in ["x", "y", "z"]:
        _r(f"bc_{d}%beg", INT, {"bc"})
        _r(f"bc_{d}%end", INT, {"bc"})

    # Relativity
    _r("relativity", LOG, {"relativity"})

    # Other (no specific feature tag)
    for n in [
        "model_eqns",
        "num_fluids",
        "thermal",
        "relax_model",
        "igr_order",
        "num_bc_patches",
        "num_patches",
        "perturb_flow_fluid",
        "perturb_sph_fluid",
        "dist_type",
        "mixlayer_perturb_nk",
        "elliptic_smoothing_iters",
        "n_start_old",
        "n_start",
        "fd_order",
        "num_igr_iters",
        "num_igr_warm_start_iters",
        "igr_iter_solver",
        "nv_uvm_igr_temps_on_gpu",
        "flux_lim",
    ]:
        _r(n, INT)
    _r("poly_sigma", REAL, math=r"\f$\sigma_\text{poly}\f$")
    _r("palpha_eps", REAL, math=r"\f$\varepsilon_\alpha\f$")
    _r("ptgalpha_eps", REAL, math=r"\f$\varepsilon_\alpha\f$")
    _r("pi_fac", REAL, math=r"\f$\pi\text{-factor}\f$")
    for n in [
        "mixlayer_vel_coef",
        "mixlayer_perturb_k0",
        "perturb_flow_mag",
        "sigR",
        "sigV",
        "rhoRV",
        "tau_star",
        "cont_damage_s",
        "alpha_bar",
        "alf_factor",
        "ic_eps",
        "ic_beta",
    ]:
        _r(n, REAL)
    for n in [
        "mpp_lim",
        "relax",
        "adv_n",
        "cont_damage",
        "igr",
        "down_sample",
        "old_grid",
        "old_ic",
        "mixlayer_vel_profile",
        "mixlayer_perturb",
        "perturb_flow",
        "perturb_sph",
        "elliptic_smoothing",
        "simplex_perturb",
        "alt_soundspeed",
        "mixture_err",
        "rdma_mpi",
        "igr_pres_lim",
        "nv_uvm_out_of_core",
        "nv_uvm_pref_gpu",
    ]:
        _r(n, LOG)
    _r("int_comp", INT)
    _r("case_dir", STR, str_len="path_len")
    _r("file_extension", STR, str_len="path_len")
    _r("files_dir", STR, str_len="path_len")

    # Body force
    for d in ["x", "y", "z"]:
        _r(f"g_{d}", REAL, math=r"\f$g_" + d + r"\f$")
        _r(f"k_{d}", REAL, math=r"\f$k_" + d + r"\f$")
        _r(f"w_{d}", REAL, math=r"\f$\omega_" + d + r"\f$")
        _r(f"p_{d}", REAL, math=r"\f$\phi_" + d + r"\f$")
        _r(f"bf_{d}", LOG)

    # Interfacial flow inputs
    _r("normMag", REAL)
    _r("p0_ic", REAL)
    _r("g0_ic", REAL)
    _r("normFac", REAL)
    _r("interface_file", STR)

    # Body force with spatial support (Wei & Freund, JFM 2005)
    _r("bf_spatial_support", LOG)
    for a in ["amp", "x_centroid", "y_centroid", "conv_vel", "sigma"]:
        _r(f"spatial_bf%{a}", REAL)
    for j in range(1, 9):
        _r(f"spatial_bf%freq({j})", REAL)
        _r(f"spatial_bf%phase({j})", REAL)

    # Synthetic turbulence forcing
    _r("synthetic_turbulence", LOG, {"synthetic_turbulence"})
    _r("synth_seed", INT, {"synthetic_turbulence"})
    _r("synth_n_shells", INT, {"synthetic_turbulence"})
    _r("num_turbulent_sources", INT, {"synthetic_turbulence"})
    _r("synth_U_inf", REAL, {"synthetic_turbulence"}, math=r"\f$U_\infty\f$")
    NSS = _fc("num_synth_shells_max", 50)
    NTS = _fc("num_turb_sources_max", 10)
    for i in range(1, NSS + 1):
        _r(f"synth_n_waves_per_shell({i})", INT, {"synthetic_turbulence"})
        _r(f"synth_k_shell({i})", REAL, {"synthetic_turbulence"}, math=r"\f$k_s\f$")
        _r(f"synth_amp_shell({i})", REAL, {"synthetic_turbulence"}, math=r"\f$A_s\f$")
    for i in range(1, NTS + 1):
        for d in range(1, 4):
            _r(f"turb_pos({i},{d})", REAL, {"synthetic_turbulence"})
            _r(f"synth_L({i},{d})", REAL, {"synthetic_turbulence"})

    # INDEXED PARAMETERS

    # patch_icpp (10 patches)
    for i in range(1, NP + 1):
        px = f"patch_icpp({i})%"
        for a in ["geometry", "smooth_patch_id", "hcid", "model_id"]:
            _r(f"{px}{a}", INT)
        for a in ["smoothen", "alter_patch"] if i >= 2 else ["smoothen"]:
            _r(f"{px}{a}", LOG)
        for a, sym in [("rho", r"\f$\rho\f$"), ("gamma", r"\f$\gamma\f$"), ("pi_inf", r"\f$\pi_\infty\f$"), ("cv", r"\f$c_v\f$"), ("qv", r"\f$q_v\f$"), ("qvp", r"\f$q'_v\f$")]:
            _r(f"{px}{a}", REAL, math=sym)
        for a in ["radius", "radii", "epsilon", "beta", "normal", "alpha_rho", "non_axis_sym", "smooth_coeff", "vel", "alpha"]:
            _r(f"{px}{a}", REAL)
        # Bubble fields
        for a in ["r0", "v0", "p0", "m0"]:
            _r(f"{px}{a}", REAL, {"bubbles"})
        for j in range(2, 10):
            _r(f"{px}a({j})", REAL)
        _r(f"{px}pres", A_REAL, math=r"\f$p\f$")
        _r(f"{px}cf_val", A_REAL)
        _r(f"{px}rxn_val", A_REAL)
        # MHD fields
        for a, sym in [("Bx", r"\f$B_x\f$"), ("By", r"\f$B_y\f$"), ("Bz", r"\f$B_z\f$")]:
            _r(f"{px}{a}", A_REAL, {"mhd"}, math=sym)
        # Chemistry species
        for j in range(1, 101):
            _r(f"{px}Y({j})", A_REAL, {"chemistry"})
        for d in ["x", "y", "z"]:
            _r(f"{px}{d}_centroid", REAL)
            _r(f"{px}length_{d}", REAL)
        for j in range(1, 4):
            _r(f"{px}radii({j})", REAL)
            _r(f"{px}normal({j})", REAL)
            _r(f"{px}vel({j})", A_REAL, math=r"\f$u_" + str(j) + r"\f$")
        for f in range(1, NF + 1):
            _r(f"{px}alpha({f})", A_REAL, math=r"\f$\alpha_" + str(f) + r"\f$")
            _r(f"{px}alpha_rho({f})", A_REAL, math=r"\f$\alpha \rho\f$")
        # Hypoelastic stress tensor
        for j in range(1, 7):
            _r(f"{px}tau_e({j})", A_REAL, {"hypoelasticity"}, math=r"\f$\tau_e\f$")
        if i >= 2:
            for j in range(1, i):
                _r(f"{px}alter_patch({j})", LOG)
        # 2D modal (geometry 13): Fourier modes and options
        for j in range(1, 11):
            _r(f"{px}fourier_cos({j})", REAL)
            _r(f"{px}fourier_sin({j})", REAL)
        _r(f"{px}modal_clip_r_to_min", LOG)
        _r(f"{px}modal_r_min", REAL)
        _r(f"{px}modal_use_exp_form", LOG)
        # 3D spherical harmonic (geometry 14): coeffs (l, m), l=0..5, m=-l..l
        for ll in range(0, 6):
            for mm in range(-ll, ll + 1):
                _r(f"{px}sph_har_coeff({ll},{mm})", REAL)

    # Derived from EOS_FAMILIES (eos_families.py), which must match the eos_* constants in src/common/m_constants.fpp.
    _EOS_NAMES = {f.suffix: f.value for f in EOS_FAMILIES}
    _EOS_VALUE_LABELS = {f.value: f.label for f in EOS_FAMILIES}
    _EOS_CHOICES = sorted(f.value for f in EOS_FAMILIES)

    # fluid_pp (10 fluids)
    # Members present in physical_parameters: gamma, pi_inf, Re, cv, qv, qvp, G.
    # mul0/ss/pv/gamma_v/M_v/mu_v/k_v/cp_v/D_v were removed from the Fortran type
    # by upstream #1085/#1093 — they must NOT be registered (namelist read would crash).
    for f in range(1, NF + 1):
        px = f"fluid_pp({f})%"
        CONSTRAINTS[f"fluid_pp({f})%eos"] = {"choices": _EOS_CHOICES, "value_labels": _EOS_VALUE_LABELS, "names": _EOS_NAMES, "fortran_prefix": "eos"}
        for a, sym in [("gamma", r"\f$\gamma_k\f$"), ("pi_inf", r"\f$\pi_{\infty,k}\f$"), ("cv", r"\f$c_{v,k}\f$"), ("qv", r"\f$q_{v,k}\f$"), ("qvp", r"\f$q'_{v,k}\f$")]:
            _r(f"{px}{a}", REAL, math=sym)
        _r(f"{px}eos", INT, math=r"\f$\mathrm{EOS}_k\f$")
        for fam in EOS_FAMILIES:
            if fam.prefix is None:
                continue
            for suffix, sym in fam.required + fam.optional:
                _r(f"{px}{fam.prefix}_{suffix}", REAL, math=sym)
        _r(f"{px}jwl_Q", REAL, math=r"\f$Q_k\f$")
        _r(f"{px}G", REAL, {"hypoelasticity"}, math=r"\f$G_k\f$")
        _r(f"{px}Re(1)", REAL, {"viscosity"}, math=r"\f$\mathrm{Re}_k\f$ (shear)")
        _r(f"{px}Re(2)", REAL, {"viscosity"}, math=r"\f$\mathrm{Re}_k\f$ (bulk)")
        _r(f"{px}k_therm", REAL, {"heat_conduction"}, math=r"\f$k_k\f$")
        _r(f"{px}non_newtonian", LOG, {"viscosity"}, math=r"\mathrm{non\text{-}Newtonian}_k")
        _r(f"{px}K", REAL, {"viscosity"}, math=r"K_k")
        _r(f"{px}nn", REAL, {"viscosity"}, math=r"n_k")
        _r(f"{px}tau0", REAL, {"viscosity"}, math=r"\tau_{0,k}")
        _r(f"{px}hb_m", REAL, {"viscosity"}, math=r"m_k")
        _r(f"{px}mu_min", REAL, {"viscosity"}, math=r"\mu_{\min,k}")
        _r(f"{px}mu_max", REAL, {"viscosity"}, math=r"\mu_{\max,k}")
        _r(f"{px}mu_bulk", REAL, {"viscosity"}, math=r"\mu_{\mathrm{bulk},k}")

    # bub_pp (bubble properties)
    for a, sym in [
        ("R0ref", r"\f$R_0\f$"),
        ("p0ref", r"\f$p_0\f$"),
        ("rho0ref", r"\f$\rho_l\f$"),
        ("T0ref", r"\f$T_0\f$"),
        ("ss", r"\f$\sigma\f$"),
        ("pv", r"\f$p_v\f$"),
        ("vd", r"\f$D\f$"),
        ("mu_l", r"\f$\mu_l\f$"),
        ("mu_v", r"\f$\mu_v\f$"),
        ("mu_g", r"\f$\mu_g\f$"),
        ("gam_v", r"\f$\gamma_v\f$"),
        ("gam_g", r"\f$\gamma_g\f$"),
        ("M_v", r"\f$M_v\f$"),
        ("M_g", r"\f$M_g\f$"),
        ("k_v", r"\f$k_v\f$"),
        ("k_g", r"\f$k_g\f$"),
        ("cp_v", r"\f$c_{p,v}\f$"),
        ("cp_g", r"\f$c_{p,g}\f$"),
        ("R_v", r"\f$R_v\f$"),
        ("R_g", r"\f$R_g\f$"),
    ]:
        _r(f"bub_pp%{a}", REAL, {"bubbles"}, math=sym)

    # patch_ib (immersed boundaries) — registered as indexed family for O(1) lookup.
    # max_index=NIB enforces the namelist limit (num_ib_patches_max_namelist); particle beds can
    # grow patch_ib beyond this at runtime, but those entries are never in the namelist.
    _ib_tags = {"ib"}
    _ib_attrs: Dict[str, tuple] = {}
    for a in ["geometry", "moving_ibm", "airfoil_id", "model_id", "inj_species"]:
        _ib_attrs[a] = (INT, _ib_tags)
    for a, pt in [("radius", REAL), ("slip", LOG), ("mass", REAL), ("v_blow", REAL), ("burn_rate_exp", REAL), ("burn_rate_pref", REAL)]:
        _ib_attrs[a] = (pt, _ib_tags)
    for j in range(1, 4):
        _ib_attrs[f"angles({j})"] = (REAL, _ib_tags)
    for d in ["x", "y", "z"]:
        _ib_attrs[f"{d}_centroid"] = (REAL, _ib_tags)
        _ib_attrs[f"length_{d}"] = (REAL, _ib_tags)
    for j in range(1, 4):
        _ib_attrs[f"vel({j})"] = (A_REAL, _ib_tags)
        _ib_attrs[f"angular_vel({j})"] = (A_REAL, _ib_tags)
    # prescribed kinematics, evaluated at run time so one binary serves every parameter value
    _ib_attrs["kin_model"] = (INT, _ib_tags)
    for j in range(1, 4):
        _ib_attrs[f"kin_hinge({j})"] = (REAL, _ib_tags)
        _ib_attrs[f"kin_offset({j})"] = (REAL, _ib_tags)
    for a in ["kin_phi0", "kin_theta0", "kin_theta_mean", "kin_freq", "kin_phase", "kin_t0", "kin_ramp", "kin_pitch_rate", "kin_smooth"]:
        _ib_attrs[a] = (REAL, _ib_tags)
    REGISTRY.register_family(
        IndexedFamily(
            base_name="patch_ib",
            attrs=_ib_attrs,
            tags=_ib_tags,
            max_index=NIB,
        )
    )

    # ib_airfoil — NACA 4-digit airfoil parameters, referenced by patch_ib(i)%airfoil_id
    _af_tags = {"ib"}
    _af_attrs: Dict[str, tuple] = {}
    for a in ["c", "p", "t", "m"]:
        _af_attrs[a] = (REAL, _af_tags)
    REGISTRY.register_family(
        IndexedFamily(
            base_name="ib_airfoil",
            attrs=_af_attrs,
            tags=_af_tags,
            max_index=NAF,
        )
    )

    # stl_models — STL/OBJ model parameters, referenced by patch_ib(i)%model_id
    _sm_attrs: Dict[str, tuple] = {}
    _sm_attrs["model_filepath"] = (STR, set())
    _sm_attrs["model_threshold"] = (REAL, set())
    for t in ["translate", "scale"]:
        for j in range(1, 4):
            _sm_attrs[f"model_{t}({j})"] = (REAL, set())
    REGISTRY.register_family(
        IndexedFamily(
            base_name="stl_models",
            attrs=_sm_attrs,
            tags=set(),
            max_index=NSM,
        )
    )

    # particle_cloud — compact bed specification that expands into individual patch_ib spheres/circles at startup
    _pb_tags = {"ib"}
    _pb_attrs: Dict[str, tuple] = {}
    for _d in ["x", "y", "z"]:
        _pb_attrs[f"{_d}_centroid"] = (REAL, _pb_tags)
        _pb_attrs[f"length_{_d}"] = (REAL, _pb_tags)
    _pb_attrs["num_particles"] = (INT, _pb_tags)
    _pb_attrs["radius"] = (REAL, _pb_tags)
    _pb_attrs["mass"] = (REAL, _pb_tags)
    _pb_attrs["min_spacing"] = (REAL, _pb_tags)
    _pb_attrs["shell_inner_radius"] = (REAL, _pb_tags)
    _pb_attrs["shell_outer_radius"] = (REAL, _pb_tags)
    _pb_attrs["moving_ibm"] = (INT, _pb_tags)
    _pb_attrs["seed"] = (INT, _pb_tags)
    _pb_attrs["cloud_geometry"] = (INT, _pb_tags)
    _pb_attrs["shell_axis"] = (INT, _pb_tags)
    _pb_attrs["packing_method"] = (INT, _pb_tags)
    _pb_attrs["periodic"] = (INT, _pb_tags)
    REGISTRY.register_family(
        IndexedFamily(
            base_name="particle_cloud",
            attrs=_pb_attrs,
            tags=_pb_tags,
            max_index=NPB,
        )
    )

    # acoustic sources (4 sources)
    for i in range(1, NA + 1):
        px = f"acoustic({i})%"
        for a in ["pulse", "support", "num_elements", "element_on", "bb_num_freq"]:
            _r(f"{px}{a}", INT, {"acoustic"})
        _r(f"{px}dipole", LOG, {"acoustic"})
        for a in [
            "mag",
            "length",
            "height",
            "wavelength",
            "frequency",
            "gauss_sigma_dist",
            "gauss_sigma_time",
            "npulse",
            "dir",
            "delay",
            "foc_length",
            "aperture",
            "element_spacing_angle",
            "element_polygon_ratio",
            "rotate_angle",
            "bb_bandwidth",
            "bb_lowest_freq",
        ]:
            _r(f"{px}{a}", REAL, {"acoustic"})
        for j in range(1, 4):
            _r(f"{px}loc({j})", REAL, {"acoustic"})

    # probes (10 probes)
    for i in range(1, NPR + 1):
        for d in ["x", "y", "z"]:
            _r(f"probe({i})%{d}", REAL, {"probes"})

    # Extended BC
    for d in ["x", "y", "z"]:
        px = f"bc_{d}%"
        for a in ["vb1", "vb2", "vb3", "ve1", "ve2", "ve3", "pres_in", "pres_out", "vel_in_ramp", "vel_in_t0", "vel_in_frac0"]:
            _r(f"{px}{a}", REAL, {"bc"})
        for a in ["grcbc_in", "grcbc_out", "grcbc_vel_out"]:
            _r(f"{px}{a}", LOG, {"bc"})
        for f in range(1, NF + 1):
            _r(f"{px}alpha_rho_in({f})", REAL, {"bc"})
            _r(f"{px}alpha_in({f})", REAL, {"bc"})
        for j in range(1, 4):
            _r(f"{px}vel_in({j})", REAL, {"bc"})
            _r(f"{px}vel_out({j})", REAL, {"bc"})

        for a in ["Twall_in", "Twall_out"]:
            _r(f"{px}{a}", REAL, {"bc"})
        for a in ["isothermal_in", "isothermal_out"]:
            _r(f"{px}{a}", LOG, {"bc"})

    # patch_bc (10 BC patches)
    for i in range(1, NB + 1):
        px = f"patch_bc({i})%"
        for a in ["geometry", "type", "dir", "loc"]:
            _r(f"{px}{a}", INT, {"bc"})
        for j in range(1, 4):
            _r(f"{px}centroid({j})", REAL, {"bc"})
            _r(f"{px}length({j})", REAL, {"bc"})
        _r(f"{px}radius", REAL, {"bc"})

    # simplex_params
    for f in range(1, NF + 1):
        _r(f"simplex_params%perturb_dens({f})", LOG)
        _r(f"simplex_params%perturb_dens_freq({f})", REAL)
        _r(f"simplex_params%perturb_dens_scale({f})", REAL)
        for j in range(1, 4):
            _r(f"simplex_params%perturb_dens_offset({f}, {j})", REAL)
    for d in range(1, 4):
        _r(f"simplex_params%perturb_vel({d})", LOG)
        _r(f"simplex_params%perturb_vel_freq({d})", REAL)
        _r(f"simplex_params%perturb_vel_scale({d})", REAL)
        for j in range(1, 4):
            _r(f"simplex_params%perturb_vel_offset({d},{j})", REAL)

    # lag_params (Lagrangian bubbles)
    # Members present in bubbles_lagrange_parameters: solver_approach, cluster_type,
    # pressure_corrector, smooth_type, heatTransfer_model, massTransfer_model,
    # write_bubbles, write_bubbles_stats, write_void_evol, pressure_force,
    # gravity_force, nBubs_glb, epsilonb, charwidth, valmaxvoid. T0/Thost/c0/rho0/x0
    # were removed from the Fortran type by upstream #1085/#1093 — they must NOT be
    # registered (namelist read would crash).
    for a in ["heatTransfer_model", "massTransfer_model", "pressure_corrector", "kahan_summation"]:
        _r(f"lag_params%{a}", LOG, {"bubbles"})
    for a in ["cluster_type", "smooth_type", "nBubs_glb"]:
        _r(f"lag_params%{a}", INT, {"bubbles"})
    _r("lag_params%charNz", INT, {"bubbles", "particles"})
    _r("lag_params%solver_approach", INT, {"bubbles", "particles"})
    for a in ["epsilonb", "valmaxvoid", "charwidth"]:
        _r(f"lag_params%{a}", REAL, {"bubbles"})
    for a in ["vel_model", "drag_model"]:
        _r(f"lag_params%{a}", INT, {"bubbles", "particles"})
    for a in ["write_bubbles", "write_bubbles_stats", "write_void_evol", "pressure_force", "gravity_force"]:
        _r(f"lag_params%{a}", LOG, {"bubbles", "particles"})
    _r("lag_params%input_path", STR, {"bubbles", "particles"})
    for a in ["nParticles_glb", "qs_drag_model", "stokes_drag", "added_mass_model", "interpolation_order", "N_collision_subcycles"]:
        _r(f"lag_params%{a}", INT, {"particles"})
    for a in ["collision_force", "subcycle_collisions", "qs_fluct_force"]:
        _r(f"lag_params%{a}", LOG, {"particles"})
    for f in range(1, NF + 1):
        _r(f"lag_params%mu_ref({f})", REAL, {"particles"})
        _r(f"lag_params%suth({f})", REAL, {"particles"})

    # chem_params
    for a in ["diffusion", "reactions", "adap_substeps"]:
        _r(f"chem_params%{a}", LOG, {"chemistry"})
    for a in ["gamma_method", "transport_model", "reaction_substeps", "reaction_substeps_max"]:
        _r(f"chem_params%{a}", INT, {"chemistry"})

    # Per-fluid output arrays
    for f in range(1, NF + 1):
        _r(f"schlieren_alpha({f})", REAL, {"output"})
        for a in ["alpha_rho_wrt", "alpha_wrt", "alpha_rho_e_wrt"]:
            _r(f"{a}({f})", LOG, {"output"})
    for j in range(1, 4):
        for a in ["mom_wrt", "vel_wrt", "flux_wrt", "omega_wrt"]:
            _r(f"{a}({j})", LOG, {"output"})

    # chem_wrt (chemistry output)
    for j in range(1, 101):
        _r(f"chem_wrt_Y({j})", LOG, {"chemistry", "output"})
    _r("chem_wrt_T", LOG, {"chemistry", "output"})

    # fluid_rho
    for f in range(1, NF + 1):
        _r(f"fluid_rho({f})", REAL)


# Load definitions when module imported and freeze registry
def _init_registry():
    """Initialize and freeze the registry. Called once at module import."""
    try:
        # Validate constraint and dependency schemas first
        # This catches typos like "choises" instead of "choices"
        _validate_all_constraints(CONSTRAINTS)
        _validate_all_dependencies(DEPENDENCIES)

        # Load all parameter definitions
        _load()

        # Freeze registry to prevent further modifications
        REGISTRY.freeze()
    except Exception as e:
        # Re-raise with context to help debugging initialization failures
        raise RuntimeError(f"Failed to initialize parameter registry: {e}\nThis is likely a bug in the parameter definitions.") from e


_init_registry()

# Namelist target mapping for Fortran codegen.
# Maps each Fortran namelist root variable to the set of MFC executables whose
# namelist it appears in. Used by fortran_gen.py to generate per-target .fpp files.
#
# When adding a new parameter:
#   1. Add to definitions.py (type, constraints, etc.) — you are here
#   2. Add the namelist root variable to NAMELIST_VARS with its target set
#   3. Re-run cmake to regenerate the .fpp files (cmake reconfigure)

NAMELIST_VARS: dict[str, set[str]] = {}

# Some common modules are compiled into every executable even though the
# corresponding user input remains meaningful in only a subset of namelists.
# Keep declaration visibility separate from input acceptance.
DECLARATION_TARGETS: dict[str, set[str]] = {}

# Maps indexed-family base names to their Fortran dimension expression.
# The generator emits `{type}, dimension({dim}) :: {name}` for each entry.
# Add here whenever a new array param needs no manual Fortran declaration.
FORTRAN_ARRAY_DIMS: dict[str, str] = {
    "fluid_rho": "num_fluids_max",
    "alpha_rho_wrt": "num_fluids_max",
    "alpha_rho_e_wrt": "num_fluids_max",
    "alpha_wrt": "num_fluids_max",
    "schlieren_alpha": "num_fluids_max",
    "chem_wrt_Y": "num_species",
    "flux_wrt": "3",
    "mom_wrt": "3",
    "omega_wrt": "3",
    "vel_wrt": "3",
    "lso_a_x": "5, 60",
    "lso_a_y": "5, 60",
    "lso_a_z": "5, 60",
    "lso2_a_x": "5, 60",
    "lso2_a_y": "5, 60",
    "lso2_a_z": "5, 60",
    "lso_pp_a_x": "5, 60",
    "lso_pp_a_y": "5, 60",
    "lso_pp_a_z": "5, 60",
    "lso_pp2_a_x": "5, 60",
    "lso_pp2_a_y": "5, 60",
    "lso_pp2_a_z": "5, 60",
}

# Derived-type namelist variables whose Fortran declarations come from generated_decls.fpp.
# Maps variable name -> (fortran_type, dimension_expr_or_None, sim_gpu_declare, doxygen_desc_or_None).
# sim_gpu_declare=True means the generator emits this variable's $:GPU_DECLARE line
# alongside its declaration; the variable must then be REMOVED from any grouped
# GPU_DECLARE list in src/simulation/m_global_parameters.fpp (multiple declare
# directives accumulate identically in OpenACC and OpenMP).
TYPED_DECLS: dict[str, tuple] = {
    "fluid_pp": ("type(physical_parameters)", "num_fluids_max", False, "Per-fluid stiffened-gas EOS parameters, Reynolds numbers, and shear modulus"),
    "bub_pp": ("type(subgrid_bubble_physical_parameters)", None, False, "Subgrid bubble physical parameters"),
    "particle_pp": ("type(subgrid_particle_physical_parameters)", None, False, "Subgrid solid-particle physical parameters"),
    "patch_icpp": ("type(ic_patch_parameters)", "num_patches_max", False, "IC patch parameters"),
    "patch_bc": ("type(bc_patch_parameters)", "num_bc_patches_max", False, "Boundary condition patch parameters"),
    "patch_ib": ("type(ib_patch_parameters)", "num_ib_patches_max_namelist", True, "Immersed boundary patch parameters"),
    "ib_airfoil": ("type(ib_airfoil_parameters)", "num_ib_airfoils_max", True, "Per-airfoil NACA user inputs"),
    "stl_models": ("type(ib_stl_parameters)", "num_stl_models_max", True, "Per-STL model parameters"),
    "probe": ("type(vec3_dt)", "num_probes_max", False, None),
    "acoustic": ("type(acoustic_parameters)", "num_probes_max", True, "Acoustic source parameters"),
    "chem_params": ("type(chemistry_parameters)", None, True, None),
    "rburn": ("type(reactive_burn_parameters)", None, True, "Condensed-phase reactive-burn (programmed detonation) parameters"),
    "lag_params": ("type(bubbles_lagrange_parameters)", None, True, "Lagrange bubbles' parameters"),
    "particle_cloud": ("type(particle_cloud_parameters)", "num_particle_clouds_max", False, "Particle bed specifications"),
    "simplex_params": ("type(simplex_noise_params)", None, False, None),
    "spatial_bf": ("type(spbf_parameters)", None, True, "Parameters for spatially supported body force (Wei & Freund, JFM 2005)"),
}


def _nv(targets: set, *names: str) -> None:
    for n in names:
        NAMELIST_VARS[n] = set(targets)


def _decl(targets: set, *names: str) -> None:
    for n in names:
        DECLARATION_TARGETS[n] = set(targets)


_ALL, _PRE_SIM, _SIM_POST = {"pre", "sim", "post"}, {"pre", "sim"}, {"sim", "post"}
_PRE_POST = {"pre", "post"}
_SIM = {"sim"}
_PRE = {"pre"}
_POST = {"post"}

_decl(
    _ALL,
    "avg_state",
    "alt_soundspeed",
    "mixture_err",
    "sigR",
    "viscous",
    "riemann_solver",
    "stl_models",
    "num_stl_models",
    "palpha_eps",
    "ptgalpha_eps",
)

_nv(
    _ALL,
    "m",
    "n",
    "p",
    "cyl_coord",
    "bc_x",
    "bc_y",
    "bc_z",
    "num_bc_patches",
    "case_dir",
    "t_step_start",
    "cfl_adap_dt",
    "cfl_const_dt",
    "n_start",
    "model_eqns",
    "mpp_lim",
    "relax",
    "relax_model",
    "fluid_pp",
    "bub_pp",
    "bubbles_euler",
    "bubbles_lagrange",
    "R0ref",
    "polytropic",
    "thermal",
    "Ca",
    "Web",
    "Re_inv",
    "polydisperse",
    "poly_sigma",
    "qbmm",
    "sigma",
    "adv_n",
    "hypoelasticity",
    "surface_tension",
    "jwl_afterburn",
    "jwl_reactive",
    "relativity",
    "ib",
    "num_ibs",
    "cont_damage",
    "hyper_cleaning",
    "Bx0",
    "precision",
    "parallel_io",
    "file_per_process",
    "fft_wrt",
    "down_sample",
)
_nv(
    _SIM_POST,
    "t_step_stop",
    "t_step_save",
    "t_stop",
    "t_save",
    "cfl_target",
    "prim_vars_wrt",
    "fd_order",
    "ib_state_wrt",
    "ib_force_wrt",
    "ib_force_stride",
    "avg_state",
    "alt_soundspeed",
    "mixture_err",
)
_nv(
    _ALL,
    "num_particle_clouds",
    "particle_cloud",
)
_nv(
    _SIM_POST,
    "lso_filter",
    "lso_filter_wrt",
    "filter_sigma",
    "lso_down_sample_factor",
    "lso_stat_wrt",
    "lso_R_gas",
    "lso_mu",
    "lso_pp_filter",
    "lso_closure_wrt",
    "lso_n_passes_x",
    "lso_n_passes_y",
    "lso_n_passes_z",
    "lso_a_x",
    "lso_a_y",
    "lso_a_z",
    "lso2_n_passes_x",
    "lso2_n_passes_y",
    "lso2_n_passes_z",
    "lso2_a_x",
    "lso2_a_y",
    "lso2_a_z",
    "lso_pp_n_passes_x",
    "lso_pp_n_passes_y",
    "lso_pp_n_passes_z",
    "lso_pp_a_x",
    "lso_pp_a_y",
    "lso_pp_a_z",
    "lso_pp2_n_passes_x",
    "lso_pp2_n_passes_y",
    "lso_pp2_n_passes_z",
    "lso_pp2_a_x",
    "lso_pp2_a_y",
    "lso_pp2_a_z",
)
_decl(
    _POST,
    "lso_filter",
    "lso_filter_wrt",
    "filter_sigma",
    "lso_down_sample_factor",
    "lso_stat_wrt",
    "lso_R_gas",
    "lso_mu",
    "lso_pp_filter",
    "lso_closure_wrt",
    "lso_n_passes_x",
    "lso_n_passes_y",
    "lso_n_passes_z",
    "lso_a_x",
    "lso_a_y",
    "lso_a_z",
    "lso2_n_passes_x",
    "lso2_n_passes_y",
    "lso2_n_passes_z",
    "lso2_a_x",
    "lso2_a_y",
    "lso2_a_z",
    "lso_pp_n_passes_x",
    "lso_pp_n_passes_y",
    "lso_pp_n_passes_z",
    "lso_pp_a_x",
    "lso_pp_a_y",
    "lso_pp_a_z",
    "lso_pp2_n_passes_x",
    "lso_pp2_n_passes_y",
    "lso_pp2_n_passes_z",
    "lso_pp2_a_x",
    "lso_pp2_a_y",
    "lso_pp2_a_z",
)
_nv(
    _PRE_SIM,
    "palpha_eps",
    "ptgalpha_eps",
    "t_step_old",
    "patch_ib",
    "pi_fac",
)
_nv(_PRE_POST, "num_fluids", "weno_order", "recon_type", "muscl_order", "mhd", "nb", "igr", "igr_order", "sigR")
_nv(_ALL, "reactive_burn", "rburn")
_nv(_PRE_SIM, "ib_airfoil")
_nv(_PRE_SIM, "stl_models", "num_stl_models")
_nv(
    _SIM,
    "dt",
    "t_step_print",
    "time_stepper",
    "adap_dt",
    "adap_dt_tol",
    "adap_dt_max_iters",
    "weno_eps",
    "teno_CT",
    "wenoz_q",
    "mp_weno",
    "weno_avg",
    "weno_Re_flux",
    "null_weights",
    "muscl_eps",
    "int_comp",
    "ic_eps",
    "ic_beta",
    "riemann_hypo_ADC",
    "ADC_kappa",
    "hll_u_interface",
    "hypo_hll_interface_rhs",
    "wave_speeds",
    "low_Mach",
    "hyper_cleaning_speed",
    "hyper_cleaning_tau",
    "run_time_info",
    "bubble_model",
    "lag_params",
    "particle_pp",
    "particles_lagrange",
    "probe_wrt",
    "num_probes",
    "jwl_ab_model",
    "jwl_q_ab",
    "jwl_ab_tau",
    "jwl_ab_A",
    "jwl_ab_theta",
    "jwl_ab_n",
    "prog_burn",
    "pb_D_cj",
    "pb_width",
    "pb_x_det",
    "pb_y_det",
    "pb_z_det",
    "pb_t_det",
    "jwl_G",
    "jwl_b_exp",
    "probe",
    "acoustic_source",
    "num_source",
    "acoustic",
    "chem_params",
    "bf_x",
    "bf_y",
    "bf_z",
    "bf_spatial_support",
    "spatial_bf",
    "k_x",
    "k_y",
    "k_z",
    "w_x",
    "w_y",
    "w_z",
    "p_x",
    "p_y",
    "p_z",
    "g_x",
    "g_y",
    "g_z",
    "synthetic_turbulence",
    "synth_seed",
    "synth_n_shells",
    "num_turbulent_sources",
    "synth_U_inf",
    "synth_n_waves_per_shell",
    "synth_k_shell",
    "synth_amp_shell",
    "turb_pos",
    "synth_L",
    "collision_model",
    "collision_temporal_resolution",
    "ramp_ratio",
    "coefficient_of_restitution",
    "collision_time",
    "ib_coefficient_of_friction",
    "ib_neighborhood_radius",
    "many_ib_patch_parallelism",
    "tau_star",
    "cont_damage_s",
    "alpha_bar",
    "rdma_mpi",
    "alf_factor",
    "num_igr_iters",
    "num_igr_warm_start_iters",
    "igr_iter_solver",
    "igr_pres_lim",
    "nv_uvm_out_of_core",
    "nv_uvm_igr_temps_on_gpu",
    "nv_uvm_pref_gpu",
    "riemann_solver",
)
_nv(
    _PRE,
    "x_domain",
    "y_domain",
    "z_domain",
    "x_a",
    "y_a",
    "z_a",
    "x_b",
    "y_b",
    "z_b",
    "stretch_x",
    "stretch_y",
    "stretch_z",
    "a_x",
    "a_y",
    "a_z",
    "loops_x",
    "loops_y",
    "loops_z",
    "n_start_old",
    "num_patches",
    "patch_icpp",
    "patch_bc",
    "sigV",
    "dist_type",
    "rhoRV",
    "old_grid",
    "old_ic",
    "perturb_flow",
    "perturb_flow_fluid",
    "perturb_flow_mag",
    "perturb_sph",
    "perturb_sph_fluid",
    "fluid_rho",
    "mixlayer_vel_profile",
    "mixlayer_vel_coef",
    "mixlayer_perturb",
    "mixlayer_perturb_nk",
    "mixlayer_perturb_k0",
    "elliptic_smoothing",
    "elliptic_smoothing_iters",
    "simplex_perturb",
    "simplex_params",
    "files_dir",
    "file_extension",
    "viscous",
)
_nv(
    _POST,
    "x_output",
    "y_output",
    "z_output",
    "format",
    "output_partial_domain",
    "sim_data",
    "G",
    "flux_lim",
    "cons_vars_wrt",
    "rho_wrt",
    "E_wrt",
    "pres_wrt",
    "c_wrt",
    "T_wrt",
    "gamma_wrt",
    "heat_ratio_wrt",
    "pi_inf_wrt",
    "pres_inf_wrt",
    "omega_wrt",
    "qm_wrt",
    "liutex_wrt",
    "schlieren_wrt",
    "schlieren_alpha",
    "alpha_rho_wrt",
    "mom_wrt",
    "vel_wrt",
    "flux_wrt",
    "alpha_wrt",
    "cf_wrt",
    "jwl_wrt",
    "chem_wrt_T",
    "chem_wrt_Y",
    "alpha_rho_e_wrt",
    "lag_header",
    "lag_txt_wrt",
    "lag_db_wrt",
    "lag_id_wrt",
    "lag_pos_wrt",
    "lag_pos_prev_wrt",
    "lag_vel_wrt",
    "lag_rad_wrt",
    "lag_rvel_wrt",
    "lag_r0_wrt",
    "lag_rmax_wrt",
    "lag_rmin_wrt",
    "lag_dphidt_wrt",
    "lag_pres_wrt",
    "lag_mv_wrt",
    "lag_mg_wrt",
    "lag_betaT_wrt",
    "lag_betaC_wrt",
)

# Case-optimization params appear in the sim namelist under #:if not MFC_CASE_OPTIMIZATION.
for _v in CASE_OPT_PARAMS:
    NAMELIST_VARS.setdefault(_v, set()).add("sim")
