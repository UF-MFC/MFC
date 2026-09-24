import copy
import dataclasses
import difflib
import json
import math
import re

import fastjsonschema

from . import common
from .analytic_expr import AnalyticExprError, fortranize_expr
from .params.eos_families import EOS_FAMILIES
from .printer import cons
from .run import case_dicts
from .state import ARG


def _suggest_similar_params(unknown_key: str, valid_keys: list, n: int = 3) -> list:
    """Find similar parameter names for typo suggestions."""
    return difflib.get_close_matches(unknown_key, valid_keys, n=n, cutoff=0.6)


QPVF_IDX_VARS = {
    "alpha_rho": "eqn_idx%cont%beg",
    "vel": "eqn_idx%mom%beg",
    "pres": "eqn_idx%E",
    "alpha": "eqn_idx%adv%beg",
    "tau_e": "eqn_idx%stress%beg",
    "Y": "eqn_idx%species%beg",
    "cf_val": "eqn_idx%c",
    "rxn_val": "eqn_idx%rxn",
    "Bx": "eqn_idx%B%beg",
    "By": "eqn_idx%B%end-1",
    "Bz": "eqn_idx%B%end",
}

# fluid_pp(i)%eos values whose coefficients depend on the local state; derived so a new
# family never leaves any_state_dependent_eos stale (see EOS_FAMILIES).
EOS_STATE_DEPENDENT_VALUES = frozenset(f.value for f in EOS_FAMILIES if f.state_dependent)

MIBM_ANALYTIC_VARS = ["vel(1)", "vel(2)", "vel(3)", "angular_vel(1)", "angular_vel(2)", "angular_vel(3)"]
# "eqn_idx%B%end - 1" not "eqn_idx%B%beg + 1" must be used because 1D does not have Bx


@dataclasses.dataclass(init=False)
class Case:
    params: dict

    def __init__(self, params: dict) -> None:
        self.params = copy.deepcopy(params)
        for key, val in self.params.items():
            self.params[key] = self.__normalize_named_value(key, val)

    @staticmethod
    def __normalize_named_value(key: str, val):
        """Convert a named enumerated value (e.g. "hllc") to its integer code."""
        from .params.errors import constraint_error
        from .params.registry import REGISTRY

        if not isinstance(val, str) or common.is_number(val):
            return val
        param_def = REGISTRY.get_param_def(key)
        if param_def is None or not param_def.constraints:
            return val
        names = param_def.constraints.get("names")
        if not names:
            return val
        if val in names:
            return names[val]
        shown = ", ".join(f"{v} ({n})" for n, v in sorted(names.items(), key=lambda kv: kv[1]))
        raise common.MFCException(constraint_error(key, "choices", shown, val))

    def get_parameters(self) -> dict:
        return self.params

    def get_cell_count(self) -> int:
        return math.prod([max(1, int(self.params.get(dir, 0))) for dir in ["m", "n", "p"]])

    def has_parameter(self, key: str) -> bool:
        return key in self.params.keys()

    def gen_json_dict_str(self) -> str:
        return json.dumps(self.params, indent=4)

    def get_inp(self, _target) -> str:
        from . import build

        target = build.get_target(_target)

        cons.print(f"Generating [magenta]{target.name}.inp[/magenta]:")
        cons.indent()

        MASTER_KEYS = case_dicts.get_input_dict_keys(target.name)

        ignored = []

        # Create Fortran-style input file content string
        dict_str = ""
        for key, val in self.params.items():
            if key in MASTER_KEYS and key not in case_dicts.IGNORE:
                if self.__is_ic_analytical(key, val) or self.__is_mib_analytical(key, val):
                    dict_str += f"{key} = 0d0\n"
                    ignored.append(key)
                    continue

                if not isinstance(val, str) or len(val) == 1:
                    dict_str += f"{key} = {val}\n"
                else:
                    dict_str += f"{key} = '{val}'\n"
            else:
                ignored.append(key)

            if key not in case_dicts.ALL:
                suggestions = _suggest_similar_params(key, list(case_dicts.ALL.keys()))
                hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
                raise common.MFCException(f"Unknown parameter '{key}'.{hint}")

        cons.print(f"[yellow]INFO:[/yellow] Forwarded {len(self.params) - len(ignored)}/{len(self.params)} parameters.")
        cons.unindent()

        # Inject the same designed weights into simulation and post_process.  Post-process
        # filtering has no IBM particle radius, so its width is defined directly by filter_sigma.
        if target.name == "simulation" and str(self.params.get("lso_filter", "F")).upper() == "T":
            dict_str += self.__get_lso_lines()
        elif target.name == "post_process" and str(self.params.get("lso_pp_filter", "F")).upper() == "T":
            dict_str += self.__get_lso_pp_lines()

        return f"&user_inputs\n{dict_str}&end/\n"

    def validate_params(self, origin_txt: str = None):
        """Validates parameters read from case file:
        1. Type checking via JSON schema
        2. Constraint validation (valid values, ranges)
        3. Dependency checking (required/recommended params)
        """
        # Type checking
        try:
            case_dicts.get_validator()(self.params)
        except fastjsonschema.JsonSchemaException as e:
            if origin_txt:
                raise common.MFCException(f"{origin_txt}: {e}")
            raise common.MFCException(f"{e}")

        # Constraint and dependency validation
        from .params.validate import validate_case

        errors, warnings = validate_case(self.params)

        # Show warnings (non-fatal)
        if warnings:
            cons.print()
            for w in warnings:
                cons.print(f"[yellow]Warning:[/yellow] {w}")

        # Raise errors (fatal)
        if errors:
            error_msg = "\n".join(f"  - {e}" for e in errors)
            if origin_txt:
                raise common.MFCException(f"{origin_txt}:\n{error_msg}")
            raise common.MFCException(f"Validation errors:\n{error_msg}")

    def __get_grid_spacing(self, down_sample_factor: int = 1):
        p = self.params
        d_p = 2.0 * float(p.get("patch_ib(1)%radius", 0.0))
        if d_p <= 0.0:
            d_p = float(p.get("filter_sigma", 0.0))
        m_cells = int(p.get("m", 0))
        n_cells = int(p.get("n", 0))
        p_cells = int(p.get("p", 0))
        if down_sample_factor > 1:
            m_cells = (m_cells + 1) // down_sample_factor - 1
            n_cells = (n_cells + 1) // down_sample_factor - 1 if n_cells > 0 else 0
            p_cells = (p_cells + 1) // down_sample_factor - 1 if p_cells > 0 else 0

        def extent(direction: str):
            return float(p.get(f"{direction}_domain%beg", 0.0)), float(p.get(f"{direction}_domain%end", 1.0))

        x_beg, x_end = extent("x")
        y_beg, y_end = extent("y")
        z_beg, z_end = extent("z")
        dx = (x_end - x_beg) / (m_cells + 1)
        dy = (y_end - y_beg) / (n_cells + 1) if n_cells > 0 else 0.0
        dz = (z_end - z_beg) / (p_cells + 1) if p_cells > 0 else 0.0
        return d_p, dx, dy, dz

    def __warn_lso_width(self, sigma: float, dx: float, dy: float, dz: float) -> None:
        for tag, spacing in (("x", dx), ("y", dy), ("z", dz)):
            if spacing > 0.0 and sigma / spacing > 40.0:
                cons.print(f"[yellow]Warning:[/yellow] LSO: sigma = {sigma / spacing:.1f} cells in {tag} is near the ~45-cell stability limit.")

    def __get_lso_lines(self) -> str:
        """Compute in-situ weights, including the optional coarse-grid second stage."""
        from .lso_filter import compute_lso_params, lso_namelist_lines

        p = self.params
        factor = int(p.get("lso_down_sample_factor", 1))
        d_p, dx, dy, dz = self.__get_grid_spacing()
        sigma = float(p.get("filter_sigma", d_p / 2.0))
        if sigma <= 0.0:
            raise common.MFCException("filter_sigma must be > 0.")
        if factor > 1 and str(p.get("lso_filter_wrt", "F")).upper() == "T":
            sigma1 = 8.0 * max(factor, 1.0) * max(d for d in (dx, dy, dz) if d > 0.0)
            if sigma > 1.05 * sigma1:
                _, cdx, cdy, cdz = self.__get_grid_spacing(down_sample_factor=factor)
                sigma2 = math.sqrt(sigma * sigma - sigma1 * sigma1)
                self.__warn_lso_width(sigma1, dx, dy, dz)
                self.__warn_lso_width(sigma2, cdx, cdy, cdz)
                lines = lso_namelist_lines(compute_lso_params(d_p, dx, dy, dz, sigma1))
                return lines + lso_namelist_lines(compute_lso_params(d_p, cdx, cdy, cdz, sigma2), prefix="lso2")
        self.__warn_lso_width(sigma, dx, dy, dz)
        return lso_namelist_lines(compute_lso_params(d_p, dx, dy, dz, sigma))

    def __get_lso_pp_lines(self) -> str:
        """Compute the post-process pass from input and target Gaussian widths."""
        from .lso_filter import PP_SPLIT_CELLS, PP_STAGE_MAX_CELLS, compute_lso_params, lso_namelist_lines

        p = self.params
        lso_filter_wrt = str(p.get("lso_filter_wrt", "F")).upper() == "T"
        factor = int(p.get("lso_down_sample_factor", 1))
        d_p, dx, dy, dz = self.__get_grid_spacing(down_sample_factor=factor if lso_filter_wrt and factor > 1 else 1)
        sigma_in = float(p.get("lso_filter_sigma_in", p.get("filter_sigma", d_p / 2.0) if lso_filter_wrt else 0.0))
        sigma_target = float(p.get("lso_filter_sigma_target", 0.0))
        if sigma_in < 0.0 or sigma_target <= sigma_in:
            raise common.MFCException("lso_filter_sigma_target must be greater than lso_filter_sigma_in >= 0.")
        sigma2 = math.sqrt(sigma_target * sigma_target - sigma_in * sigma_in)
        d_min = min(d for d in (dx, dy, dz) if d > 0.0)
        if sigma2 / d_min > PP_SPLIT_CELLS:
            sigma_stage = sigma2 / math.sqrt(2.0)
            if sigma_stage / d_min > PP_STAGE_MAX_CELLS:
                raise common.MFCException(
                    f"lso_pp_filter: sigma_2 = {sigma2 / d_min:.1f} cells exceeds what the two post_process "
                    f"cascades can carry (~{PP_STAGE_MAX_CELLS:.0f} cells each). Filter in situ with "
                    f"lso_down_sample_factor > 1 so the post_process pass runs on the coarse grid, or lower "
                    f"lso_filter_sigma_target."
                )
            cons.print(
                f"[cyan]LSO filter (post_process):[/cyan] sigma_in={sigma_in:.4g}, "
                f"sigma_target={sigma_target:.4g}, sigma2={sigma2:.4g} as two cascades of "
                f"{sigma_stage:.4g} ({sigma_stage / d_min:.1f} cells each), computing weights..."
            )
            stage = compute_lso_params(d_p, dx, dy, dz, sigma_stage)
            return lso_namelist_lines(stage, prefix="lso_pp") + lso_namelist_lines(stage, prefix="lso_pp2")
        self.__warn_lso_width(sigma2, dx, dy, dz)
        return lso_namelist_lines(compute_lso_params(d_p, dx, dy, dz, sigma2), prefix="lso_pp")

    def __get_ndims(self) -> int:
        return 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)

    def __is_ic_analytical(self, key: str, val: str) -> bool:
        """Is this initial condition analytical?
        More precisely, is this an arbitrary expression or a string representing a number?"""
        if common.is_number(val) or not isinstance(val, str):
            return False

        for array in QPVF_IDX_VARS:
            if re.match(rf"^patch_icpp\([0-9]+\)%{array}", key):
                return True

        return False

    def __is_mib_analytical(self, key: str, val: str) -> bool:
        """Is this initial condition analytical?
        More precisely, is this an arbitrary expression or a string representing a number?"""
        if common.is_number(val) or not isinstance(val, str):
            return False

        for variable in MIBM_ANALYTIC_VARS:
            if re.match(rf"^patch_ib\([0-9]+\)%{re.escape(variable)}", key):
                return True

        return False

    def __get_analytic_ic_fpp(self, print: bool) -> str:
        # generates the content of an FFP file that will hold the functions for
        # some initial condition
        DATA = {
            1: {"ptypes": [1, 15, 16], "sf_idx": "i, 0, 0"},
            2: {"ptypes": [2, 3, 4, 5, 6, 7, 13, 17, 18, 21], "sf_idx": "i, j, 0"},
            3: {"ptypes": [8, 9, 10, 11, 12, 14, 19, 21], "sf_idx": "i, j, k"},
        }[self.__get_ndims()]

        patches = {}

        # iterates over the parameters and checks if they are defined as an
        # analytical function. If so, append it to the `patches`` object
        for key, val in self.params.items():
            if not self.__is_ic_analytical(key, val):
                continue

            patch_id = re.search(r"[0-9]+", key).group(0)

            if patch_id not in patches:
                patches[patch_id] = []

            patches[patch_id].append((key, val))

        srcs = []

        # for each analytical patch that is required to be added, generate
        # the string that contains that function.
        for pid, items in patches.items():
            ptype = self.params[f"patch_icpp({pid})%geometry"]

            if ptype not in DATA["ptypes"]:
                raise common.MFCException(f"Patch #{pid} of type {ptype} cannot be analytically defined.")

            var_map = {
                "x": "x_cc(i)",
                "y": "y_cc(j)",
                "z": "z_cc(k)",
                "xc": f"patch_icpp({pid})%x_centroid",
                "yc": f"patch_icpp({pid})%y_centroid",
                "zc": f"patch_icpp({pid})%z_centroid",
                "lx": f"patch_icpp({pid})%length_x",
                "ly": f"patch_icpp({pid})%length_y",
                "lz": f"patch_icpp({pid})%length_z",
                "r": f"patch_icpp({pid})%radius",
                "eps": f"patch_icpp({pid})%epsilon",
                "beta": f"patch_icpp({pid})%beta",
                "tau_e": f"patch_icpp({pid})%tau_e",
                "radii": f"patch_icpp({pid})%radii",
            }

            lines = []
            # translate each analytic expression to Fortran and emit an assignment
            for attribute, expr in items:
                if print:
                    cons.print(f"* Codegen: {attribute} = {expr}")

                varname = re.findall(r"[a-zA-Z][a-zA-Z0-9_]*", attribute)[1]
                qpvf_idx = QPVF_IDX_VARS[varname][:]

                if len(re.findall(r"[0-9]+", attribute)) == 2:
                    idx = int(re.findall(r"[0-9]+", attribute)[1]) - 1
                    qpvf_idx = f"{qpvf_idx} + {idx}"

                lhs = f"q_prim_vf({qpvf_idx})%sf({DATA['sf_idx']})"
                try:
                    rhs = fortranize_expr(expr, var_map)
                except AnalyticExprError as err:
                    raise common.MFCException(f"{attribute}: {err}") from err

                lines.append(f"        {lhs} = {rhs}")

            # concatenates all of the analytic lines into a single string with
            # each element separated by new line characters. Then write those
            # new lines as a fully concatenated string with fortran syntax
            srcs.append(f"""\
    if (patch_id == {pid}) then
{f"{chr(10)}".join(lines)}
    end if\
""")

        content = f"""\
! This file was generated by MFC. It is only used when analytical patches are
! present in the input file. It is used to define the analytical patches with
! expressions that are evaluated at runtime from the input file.

#:def analytical()
{f"{chr(10)}{chr(10)}".join(srcs)}
#:enddef
"""
        return content

    # gets the analytic description of a moving IB's velocity and rotation rate
    def __get_analytic_mib_fpp(self, print: bool) -> str:
        # iterates over the parameters and checks if they are defined as an
        # analytical function. If so, append it to the `patches`` object
        ib_patches = {}

        for key, val in self.params.items():
            if not self.__is_mib_analytical(key, val):
                continue

            patch_id = re.search(r"[0-9]+", key).group(0)

            if patch_id not in ib_patches:
                ib_patches[patch_id] = []

            ib_patches[patch_id].append((key, val))

        srcs = []

        # for each analytical patch that is required to be added, generate
        # the string that contains that function.
        for pid, items in ib_patches.items():
            var_map = {
                "x": "x_cc(i)",
                "y": "y_cc(j)",
                "z": "z_cc(k)",
                "t": "mytime",
                "r": f"patch_ib({pid})%radius",
            }

            lines = []
            # translate each analytic expression to Fortran and emit an assignment
            for attribute, expr in items:
                if print:
                    cons.print(f"* Codegen: {attribute} = {expr}")

                lhs = attribute
                try:
                    rhs = fortranize_expr(expr, var_map)
                except AnalyticExprError as err:
                    raise common.MFCException(f"{attribute}: {err}") from err

                lines.append(f"        {lhs} = {rhs}")

            # concatenates all of the analytic lines into a single string with
            # each element separated by new line characters. Then write those
            # new lines as a fully concatenated string with fortran syntax
            srcs.append(f"""\
    if (gbl_id == {pid}) then
{f"{chr(10)}".join(lines)}
    end if\
""")

        content = f"""\
! This file was generated by MFC. It is only used when we analytically
! parameterize the velocity and rotation rate of a moving IB.

#:def mib_analytical()
gbl_id = patch_ib(i)%gbl_patch_id
{f"{chr(10)}{chr(10)}".join(srcs)}
#:enddef
"""
        return content

    def __get_sim_fpp(self, print: bool) -> str:
        if ARG("case_optimization"):
            if print:
                cons.print("Case optimization is enabled.")

            nterms = 1

            bubble_model = int(self.params.get("bubble_model", "-100"))

            if bubble_model == 2:
                nterms = 32
            elif bubble_model == 3:
                nterms = 7

            mapped_weno = 1 if self.params.get("mapped_weno", "F") == "T" else 0
            wenoz = 1 if self.params.get("wenoz", "F") == "T" else 0
            teno = 1 if self.params.get("teno", "F") == "T" else 0
            wenojs = 0 if (mapped_weno or wenoz or teno) else 1
            recon_type = self.params.get("recon_type", 1)

            # This fixes a bug on Frontier to do with allocating 0:0 arrays
            weno_order = int(self.params.get("weno_order", 0))
            if recon_type == 1:
                weno_polyn = int((weno_order - 1) / 2)
            else:
                weno_polyn = 1

            if self.params.get("igr", "F") == "T":
                weno_order = 5
                weno_polyn = 3

            if teno:
                weno_num_stencils = weno_order - 3
            else:
                weno_num_stencils = weno_polyn

            num_dims = 1 + min(int(self.params.get("n", 0)), 1) + min(int(self.params.get("p", 0)), 1)
            if self.params.get("mhd", "F") == "T":
                num_vels = 3
            else:
                num_vels = num_dims

            # Baking this lets the compiler drop the state-dependent EOS chain entirely. Left in the
            # call graph it costs registers, and so occupancy, in every kernel that can reach it.
            eos_state_dependent = EOS_STATE_DEPENDENT_VALUES
            num_fluids_case = int(self.params.get("num_fluids", 1))
            any_state_dependent_eos = 1 if any(int(self.params.get(f"fluid_pp({f})%eos", 1)) in eos_state_dependent for f in range(1, num_fluids_case + 1)) else 0

            mhd = 1 if self.params.get("mhd", "F") == "T" else 0
            relativity = 1 if self.params.get("relativity", "F") == "T" else 0
            viscous = 1 if self.params.get("viscous", "F") == "T" else 0
            igr = 1 if self.params.get("igr", "F") == "T" else 0
            igr_pres_lim = 1 if self.params.get("igr_pres_lim", "F") == "T" else 0

            # Throw error if wenoz_q is required but not set
            out = f"""\
#:set MFC_CASE_OPTIMIZATION = {ARG("case_optimization")}
#:set recon_type            = {recon_type}
#:set weno_order            = {weno_order}
#:set weno_polyn            = {weno_polyn}
#:set muscl_order           = {int(self.params.get("muscl_order", 1))}
#:set muscl_polyn           = {int(self.params.get("muscl_order", 1))}
#:set muscl_lim             = {int(self.params.get("muscl_lim", 1))}
#:set weno_num_stencils     = {weno_num_stencils}
#:set nb                    = {int(self.params.get("nb", 1))}
#:set num_dims              = {num_dims}
#:set num_vels              = {num_vels}
#:set nterms                = {nterms}
#:set num_fluids            = {int(self.params.get("num_fluids", 1))}
#:set wenojs                = {wenojs}
#:set mapped_weno           = {mapped_weno}
#:set wenoz                 = {wenoz}
#:set teno                  = {teno}
#:set wenoz_q               = {self.params.get("wenoz_q", -1)}
#:set mhd                   = {mhd}
#:set relativity            = {relativity}
#:set igr                   = {igr}
#:set igr_iter_solver       = {self.params.get("igr_iter_solver", 1)}
#:set igr_pres_lim          = {igr_pres_lim}
#:set igr_order             = {self.params.get("igr_order", 3)}
#:set viscous               = {viscous}
#:set any_state_dependent_eos = {any_state_dependent_eos}
"""

        else:
            out = ""

        out = out + self.__get_analytic_mib_fpp(print)

        # We need to also include the pre_processing includes so that common subroutines have access to the @:analytical function
        return out + f"\n{self.__get_pre_fpp(print)}"

    def __get_pre_fpp(self, print: bool) -> str:
        out = self.__get_analytic_ic_fpp(print)
        return out

    def get_fpp(self, target, print=True) -> str:
        from . import build

        def _prepend() -> str:
            return f"""\
#:set chemistry             = {self.params.get("chemistry", "F") == "T"}
"""

        def _default(_) -> str:
            return "! This file is purposefully empty."

        result = {
            "pre_process": self.__get_pre_fpp,
            "simulation": self.__get_sim_fpp,
        }.get(build.get_target(target).name, _default)(print)

        return _prepend() + result

    def __getitem__(self, key: str) -> str:
        return self.params[key]

    def __setitem__(self, key: str, val: str):
        self.params[key] = self.__normalize_named_value(key, val)
