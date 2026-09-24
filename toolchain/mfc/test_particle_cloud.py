import json
import os
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest

from mfc.case_validator import CaseConstraintError, CaseValidator


def _valid_cloud_params():
    return {
        "m": 31,
        "n": 31,
        "p": 31,
        "model_eqns": 2,
        "num_fluids": 1,
        "num_patches": 1,
        "t_step_start": 0,
        "t_step_stop": 1,
        "t_step_save": 1,
        "dt": 1.0e-6,
        "x_domain%beg": -1.0,
        "x_domain%end": 1.0,
        "y_domain%beg": -1.0,
        "y_domain%end": 1.0,
        "z_domain%beg": -1.0,
        "z_domain%end": 1.0,
        "bc_x%beg": -1,
        "bc_x%end": -1,
        "bc_y%beg": -1,
        "bc_y%end": -1,
        "bc_z%beg": -1,
        "bc_z%end": -1,
        "riemann_solver": 2,
        "wave_speeds": 1,
        "avg_state": 2,
        "fd_order": 2,
        "ib": "T",
        "num_ibs": 0,
        "num_particle_clouds": 1,
        "particle_cloud(1)%cloud_geometry": 2,
        "particle_cloud(1)%packing_method": 1,
        "particle_cloud(1)%x_centroid": 0.0,
        "particle_cloud(1)%y_centroid": 0.0,
        "particle_cloud(1)%z_centroid": 0.0,
        "particle_cloud(1)%num_particles": 8,
        "particle_cloud(1)%radius": 0.04,
        "particle_cloud(1)%mass": 1.0,
        "particle_cloud(1)%min_spacing": 0.01,
        "particle_cloud(1)%shell_inner_radius": 0.1,
        "particle_cloud(1)%shell_outer_radius": 0.3,
    }


def test_hemi_shell_validator_accepts_feasible_shell():
    CaseValidator(_valid_cloud_params()).validate("simulation")


def test_hemi_shell_validator_rejects_negative_inner_radius():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_inner_radius": -0.01}
    with pytest.raises(CaseConstraintError, match="shell_inner_radius >= 0"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_thinner_than_particle_diameter():
    params = {**_valid_cloud_params(), "particle_cloud(1)%shell_outer_radius": 0.17}
    with pytest.raises(CaseConstraintError, match="shell_outer_radius > shell_inner_radius"):
        CaseValidator(params).validate("simulation")


def test_hemi_shell_validator_rejects_shell_outside_domain():
    params = {**_valid_cloud_params(), "particle_cloud(1)%x_centroid": 0.9, "particle_cloud(1)%shell_outer_radius": 0.3}
    with pytest.raises(CaseConstraintError, match="x-extent must lie within x_domain"):
        CaseValidator(params).validate("simulation")


@pytest.mark.skipif(os.environ.get("MFC_TEST_CLOUD_OUTPUT") != "1", reason="requires prebuilt CPU MFC executables")
@pytest.mark.parametrize("ranks,file_per_process,shell_axis", [(1, False, 1), (2, False, 2), (1, True, 2), (2, True, 1)])
def test_simulation_reads_preprocessed_cloud(tmp_path, ranks, file_per_process, shell_axis):
    """Simulation and restart must retain the placement written by pre_process."""
    from .test.case import BASE_CFG

    root = Path(__file__).resolve().parents[2]
    params = {k: v for k, v in BASE_CFG.items() if not k.startswith(("patch_icpp(2)", "patch_icpp(3)"))}
    params.update(_valid_cloud_params())
    params.pop("bc_z%beg")
    params.pop("bc_z%end")
    params.update(
        {
            "m": 63,
            "n": 63,
            "p": 0,
            "parallel_io": "T",
            "file_per_process": "T" if file_per_process else "F",
            "ib_state_wrt": "T",
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.0,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": 2.0,
            "patch_icpp(1)%length_y": 2.0,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "particle_cloud(1)%num_particles": 4,
            "particle_cloud(1)%radius": 0.08,
            "particle_cloud(1)%shell_outer_radius": 0.6,
            "particle_cloud(1)%shell_axis": shell_axis,
            "particle_cloud(1)%seed": 12345,
            "particle_cloud(1)%moving_ibm": 0,
        }
    )
    case_path = tmp_path / "case.json"

    def run(target):
        case_path.write_text(json.dumps(params))
        result = subprocess.run(
            [str(root / "mfc.sh"), "run", str(case_path), "-t", target, "-n", str(ranks), "--no-build", "--no-gpu", "--no-debug", "--no-reldebug"],
            cwd=root,
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    def positions(step):
        records = {}
        if file_per_process:
            for rank in range(ranks):
                data = (tmp_path / f"restart_data/lustre_{step}/ib_state_{step}_{rank:07d}.dat").read_bytes()
                count = struct.unpack_from("=i", data)[0]
                assert len(data) == 4 + count * 164
                for idx in range(count):
                    offset = 4 + idx * 164
                    patch_id = struct.unpack_from("=i", data, offset)[0]
                    assert patch_id not in records
                    records[patch_id] = np.frombuffer(data, dtype=np.float64, count=20, offset=offset + 4)[16:20]
        else:
            data = np.fromfile(tmp_path / f"restart_data/ib_state_{step}.dat", dtype=np.float64).reshape(-1, 20)
            records = {idx + 1: row[16:20] for idx, row in enumerate(data)}
        assert set(records) == {1, 2, 3, 4}
        return np.array([records[idx] for idx in sorted(records)])

    run("pre_process")
    initial = positions(0)
    params["particle_cloud(1)%seed"] = 54321
    run("simulation")
    np.testing.assert_array_equal(positions(0), initial)
    np.testing.assert_array_equal(positions(1), initial)
    params.update({"t_step_start": 1, "t_step_stop": 2})
    run("simulation")
    np.testing.assert_array_equal(positions(2), initial)
