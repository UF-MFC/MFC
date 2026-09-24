"""Opt-in executable regression: build first, then set MFC_TEST_LSO_OUTPUT=1 for ./mfc.sh lint."""

import contextlib
import io
import json
import os
import runpy
import subprocess
from pathlib import Path

import numpy as np
import pytest


@pytest.mark.skipif(os.environ.get("MFC_TEST_LSO_OUTPUT") != "1", reason="requires prebuilt MFC executables")
@pytest.mark.parametrize(
    "stiffening,boundary,factor,mode",
    [(s, b, f, "fixed") for s in (0.0, 1.0) for b in (-1, -3) for f in (1, 2)] + [(0.0, -1, 2, mode) for mode in ("cfl_adap_dt", "cfl_const_dt")],
)
def test_viscous_closure_uses_each_cells_mask(tmp_path, stiffening, boundary, factor, mode):
    from .viz.silo_reader import read_silo_file

    root = Path(__file__).resolve().parents[2]
    with contextlib.redirect_stdout(io.StringIO()) as output:
        runpy.run_path(str(root / "examples/1D_lso_mach3_shocktube/case.py"))
    case = {k: v for k, v in json.loads(output.getvalue()).items() if not k.startswith("patch_icpp(2)")}
    case.update(
        {
            "m": 31,
            "n": 31,
            "dt": 1e-5,
            "t_step_stop": 2,
            "t_step_save": 2,
            "y_domain%beg": 0.0,
            "y_domain%end": 1.0,
            "bc_x%beg": boundary,
            "bc_x%end": boundary,
            "bc_y%beg": boundary,
            "bc_y%end": boundary,
            "num_patches": 1,
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%y_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%length_y": 1.0,
            "patch_icpp(1)%vel(1)": 0.2,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 1.0,
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "ib": "T",
            "num_ibs": 1,
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": 0.5,
            "patch_ib(1)%y_centroid": 0.5,
            "patch_ib(1)%radius": 0.15,
            "patch_ib(1)%slip": "F",
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 100.0,
            "lso_mu": 0.01,
            "filter_sigma": 0.06,
            "lso_pp_filter": "F",
            "lso_down_sample_factor": factor,
        }
    )
    if stiffening:
        case.update({"fluid_pp(1)%eos": "stiffened_gas", "fluid_pp(1)%pi_inf": 3.5 * stiffening})
    if mode == "wide":
        case["filter_sigma"] = 0.6
        case.update({"bc_x%beg": -1, "bc_x%end": -1})
    if mode.startswith("cfl_"):
        for key in ("dt", "t_step_start", "t_step_stop", "t_step_save"):
            case.pop(key)
        case.update({mode: "T", "n_start": 0, "t_stop": 2e-4, "t_save": 1e-4, "cfl_target": 0.002})
    case_path = tmp_path / "case.json"
    case_path.write_text(json.dumps(case))
    result = subprocess.run(
        [str(root / "mfc.sh"), "run", str(case_path), "-n", "1", "--no-debug", "--no-build"],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    if mode.startswith("cfl_"):
        marker_path = tmp_path / "restart_data/lustre_ib.dat"
        initial = np.fromfile(marker_path, dtype=np.int32, count=32 * 32)
        assert np.any(initial == 1) and np.any(initial == 0)
        for step in (1, 2):
            saved = np.fromfile(marker_path, dtype=np.int32, count=32 * 32, offset=32 * 32 * 8 * (5 + step))
            np.testing.assert_array_equal(saved, initial, err_msg=f"IBM save index {step}")
    assert result.returncode == 0, result.stdout + result.stderr
    fields = read_silo_file(str(tmp_path / "silo_hdf5/p0/2.silo")).variables
    stats = [fields[f"lso_stat_{i:03}"] for i in range(1, 22)]
    closures = [fields[f"lso_closure_{i:03}"] for i in range(1, 18)]
    assert all(np.isfinite(field).all() for field in stats + closures)
    cells = 32 // factor
    assert stats[0].size == cells * cells
    markers = np.fromfile(tmp_path / "restart_data/lustre_ib.dat", dtype=np.int32, count=32 * 32).reshape((32, 32), order="F")
    first = (factor - 1) // 2
    np.testing.assert_array_equal(fields["ib_markers"], markers[first::factor, first::factor])
    weight = np.fromfile(tmp_path / "restart_data/lustre_lso_mask_2.dat").reshape((cells, cells), order="F")
    np.testing.assert_allclose(weight + stats[0], 1.0, rtol=0, atol=1e-10)
    usable = (weight > 0.1) & (weight < 0.9)
    assert np.count_nonzero(usable) > 4, "test must exercise a varying gas mask"
    rho, mx, my = stats[1], stats[5], stats[6]
    expected_temperature = (fields["E"][usable] * weight[usable] - 0.5 * stats[2][usable] - stiffening * weight[usable]) / (2.5 * rho[usable])
    np.testing.assert_allclose(closures[14][usable], expected_temperature, rtol=1e-9, atol=1e-12)
    for actual, product, contraction in (
        (closures[7], stats[19], stats[14] * mx + stats[15] * my),
        (closures[8], stats[20], stats[15] * mx + stats[16] * my),
    ):
        expected = (product[usable] - contraction[usable] / weight[usable]) / rho[usable]
        assert np.max(np.abs(contraction[usable])) > 1e-8, "zero stress cannot detect bad normalization"
        np.testing.assert_allclose(actual[usable], expected, rtol=1e-9, atol=1e-12)


@pytest.mark.skipif(os.environ.get("MFC_TEST_LSO_OUTPUT") != "1", reason="requires prebuilt MFC executables")
@pytest.mark.parametrize(
    "ndim,factor,ranks,mode",
    [(1, 2, 1, "fixed"), (1, 3, 2, "fixed"), (2, 2, 2, "fixed"), (3, 2, 2, "fixed")] + [(1, 2, 2, mode) for mode in ("cfl_adap_dt", "cfl_const_dt", "unaligned", "wide")] + [(1, 1, 1, "serial")],
)
def test_lso_reader_preserves_fields_and_cell_centres(tmp_path, ndim, factor, ranks, mode):
    from .viz.silo_reader import read_silo_file

    root = Path(__file__).resolve().parents[2]
    with contextlib.redirect_stdout(io.StringIO()) as output:
        runpy.run_path(str(root / "examples/1D_lso_mach3_shocktube/case.py"))
    case = json.loads(output.getvalue())
    cells = 32 * factor
    if mode == "unaligned":
        cells = 90
    case.update(
        {
            "m": cells - 1,
            "weno_order": 3,
            "mp_weno": "F",
            "dt": 1e-8,
            "t_step_stop": 2,
            "t_step_save": 2,
            "filter_sigma": 0.04,
            "lso_pp_filter": "F",
            "lso_down_sample_factor": factor,
            "cons_vars_wrt": "T",
            "prim_vars_wrt": "F",
        }
    )
    for axis, key in zip("yz"[: ndim - 1], ("n", "p")):
        case.update({key: cells - 1, f"{axis}_domain%beg": 0.0, f"{axis}_domain%end": 1.0, f"bc_{axis}%beg": -3, f"bc_{axis}%end": -3})
        for patch in (1, 2):
            case.update({f"patch_icpp({patch})%{axis}_centroid": 0.5, f"patch_icpp({patch})%length_{axis}": 1.0, f"patch_icpp({patch})%vel({'xyz'.index(axis) + 1})": 0.03 * patch * "xyz".index(axis)})
        length = {"y": 0.765, "z": 0.615}[axis]
        case.update({f"patch_icpp(2)%{axis}_centroid": 1 - length / 2, f"patch_icpp(2)%length_{axis}": length})
    for patch, centre, length in ((1, 0.5, 1.0), (2, 0.7425, 0.515)):
        case.update({f"patch_icpp({patch})%geometry": {1: 1, 2: 3, 3: 9}[ndim], f"patch_icpp({patch})%x_centroid": centre, f"patch_icpp({patch})%length_x": length})
    case["patch_icpp(2)%alter_patch(1)"] = "T"
    if mode == "serial":
        case.update({"parallel_io": "F", "lso_stat_wrt": "F", "lso_closure_wrt": "F"})
    if mode == "wide":
        case.update({"filter_sigma": 0.3, "bc_x%beg": -1, "bc_x%end": -1})
    if mode.startswith("cfl_"):
        for key in ("dt", "t_step_start", "t_step_stop", "t_step_save"):
            case.pop(key)
        case.update({mode: "T", "n_start": 0, "t_stop": 2e-4, "t_save": 1e-4, "cfl_target": 0.05})
    case_path = tmp_path / "case.json"
    case_path.write_text(json.dumps(case))
    options = ["--debug", "-j", "4"] if ndim == 3 else ["--no-debug", "--no-build"]
    result = subprocess.run([str(root / "mfc.sh"), "run", str(case_path), "-n", str(ranks), *options], cwd=root, capture_output=True, text=True, timeout=180, check=False)
    if mode == "unaligned":
        assert result.returncode != 0
        assert "LSO downsampling requires per-rank cells" in result.stdout + result.stderr
        assert not (tmp_path / "restart_data/lustre_lso_2.dat").exists()
        return
    assert result.returncode == 0, result.stdout + result.stderr
    names = ["alpha_rho1", *[f"mom{i}" for i in range(1, ndim + 1)], "E", "alpha1"]
    coarse = cells // factor
    for step in (0, 2):
        count = cells if step == 0 else coarse
        prefix = "" if step == 0 else "lso_"
        if mode == "serial":
            raw = np.stack([np.fromfile(tmp_path / f"p_all/p0/{step}/{prefix}q_cons_vf{i}.dat", offset=4, count=cells) for i in range(1, 5)], axis=-1)
            if step > 0:
                unfiltered = np.fromfile(tmp_path / f"p_all/p0/{step}/q_cons_vf1.dat", offset=4, count=cells)
                assert np.max(np.abs(unfiltered - raw[:, 0])) > 1e-3
        else:
            raw = np.fromfile(tmp_path / f"restart_data/lustre_{prefix}{step}.dat").reshape((count,) * ndim + (len(names),), order="F")
        if mode == "wide" and step > 0:
            stats = np.fromfile(tmp_path / f"restart_data/lustre_lso_stat_{step}.dat")
            assert stats.size == coarse * 11
            assert np.isfinite(stats).all()
        if step == 0:
            for axis in range(ndim):
                positions = (np.arange(coarse) + 0.5) * factor - 0.5
                lower = np.floor(positions).astype(int)
                raw = 0.5 * (np.take(raw, lower, axis=axis) + np.take(raw, np.ceil(positions).astype(int), axis=axis))
        for rank in range(ranks):
            data = read_silo_file(str(tmp_path / f"silo_hdf5/p{rank}/{step}.silo"))
            coords = [data.x_cb, data.y_cb, data.z_cb][:ndim]
            indices, valid = [], []
            for cb in coords:
                np.testing.assert_allclose(np.diff(cb), 1 / coarse, rtol=0, atol=1e-12)
                idx = np.rint((cb[1:] + cb[:-1]) * 0.5 * coarse - 0.5).astype(int)
                inside = (idx >= 0) & (idx < coarse)
                indices.append(idx[inside])
                valid.append(np.flatnonzero(inside))
            for i, name in enumerate(names):
                actual = data.variables[name].reshape(tuple(len(cb) - 1 for cb in coords), order="F")
                np.testing.assert_allclose(actual[np.ix_(*valid)], raw[..., i][np.ix_(*indices)], rtol=1e-12, atol=1e-12, err_msg=f"step={step}, rank={rank}, field={name}")


@pytest.mark.skipif(os.environ.get("MFC_TEST_LSO_OUTPUT") != "1", reason="requires an MFC build")
def test_lso_linear_velocity_has_no_viscous_subfilter_stress(tmp_path):
    from .viz.silo_reader import read_silo_file

    root = Path(__file__).resolve().parents[2]
    with contextlib.redirect_stdout(io.StringIO()) as output:
        runpy.run_path(str(root / "examples/1D_lso_mach3_shocktube/case.py"))
    case = {k: v for k, v in json.loads(output.getvalue()).items() if not k.startswith("patch_icpp(2)")}
    case.update(
        {
            "m": 63,
            "dt": 1e-8,
            "t_step_stop": 2,
            "t_step_save": 2,
            "num_patches": 1,
            "patch_icpp(1)%x_centroid": 0.5,
            "patch_icpp(1)%length_x": 1.0,
            "patch_icpp(1)%vel(1)": "0.1*x",
            "patch_icpp(1)%alpha_rho(1)": 1.0,
            "patch_icpp(1)%pres": 1.0,
            "viscous": "T",
            "fluid_pp(1)%Re(1)": 100.0,
            "lso_mu": 0.01,
            "filter_sigma": 0.02,
            "lso_pp_filter": "F",
            "lso_down_sample_factor": 1,
        }
    )
    path = tmp_path / "case.json"
    path.write_text(json.dumps(case))
    result = subprocess.run([str(root / "mfc.sh"), "run", str(path), "-n", "1", "-j", "4", "--no-debug"], cwd=root, capture_output=True, text=True, timeout=180, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    fields = read_silo_file(str(tmp_path / "silo_hdf5/p0/2.silo")).variables
    np.testing.assert_allclose(fields["lso_stat_009"][20:44], 4 / 3 * 0.01 * 0.1, rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(fields["lso_closure_005"][20:44], 0.0, rtol=0, atol=1e-10)


@pytest.mark.skipif(os.environ.get("MFC_TEST_LSO_OUTPUT") != "1", reason="requires prebuilt MFC executables")
@pytest.mark.parametrize("ranks,conductivity,boundary", [(1, 0.0, -3), (2, 0.01, -1)])
def test_lso_temperature_products_use_saved_state(tmp_path, ranks, conductivity, boundary):
    root = Path(__file__).resolve().parents[2]
    with contextlib.redirect_stdout(io.StringIO()) as output:
        runpy.run_path(str(root / "examples/1D_lso_mach3_shocktube/case.py"))
    case = json.loads(output.getvalue())
    case.update(
        {
            "m": 63,
            "weno_order": 3,
            "mp_weno": "F",
            "dt": 1e-4,
            "t_step_stop": 4,
            "t_step_save": 1,
            "filter_sigma": 0.04,
            "lso_pp_filter": "F",
            "lso_down_sample_factor": 1,
            "fluid_pp(1)%k_therm": conductivity,
            "bc_x%beg": boundary,
            "bc_x%end": boundary,
        }
    )
    for patch, centre, length in ((1, 0.5, 1.0), (2, 0.7425, 0.515)):
        case.update({f"patch_icpp({patch})%x_centroid": centre, f"patch_icpp({patch})%length_x": length})
    case["patch_icpp(2)%alter_patch(1)"] = "T"
    path = tmp_path / "case.json"
    path.write_text(json.dumps(case))
    # Clear analytic IC code compiled by the linear-velocity regression.
    result = subprocess.run([str(root / "mfc.sh"), "run", str(path), "-n", str(ranks), "--no-debug", "-j", "4"], cwd=root, capture_output=True, text=True, timeout=180, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    initial = np.fromfile(tmp_path / "restart_data/lustre_0.dat").reshape((4, 64))
    np.testing.assert_allclose(initial[1, 0] / initial[0, 0], case["patch_icpp(1)%vel(1)"], rtol=1e-14)
    inp = dict(line.split(" = ", 1) for line in (tmp_path / "simulation.inp").read_text().splitlines() if " = " in line)
    passes = int(inp["lso_n_passes_x"])
    coefficients = np.array([float(v.replace("d", "e")) for v in inp["lso_a_x"].split()[: 5 * passes]]).reshape((-1, 5))
    for step in (1, 4):
        rho, momentum, energy, _ = np.fromfile(tmp_path / f"restart_data/lustre_{step}.dat").reshape((4, 64))
        temperature = (energy - 0.5 * momentum**2 / rho) / (case["fluid_pp(1)%gamma"] * case["lso_R_gas"] * rho)
        flux = -conductivity * (np.roll(temperature, -1) - np.roll(temperature, 1)) / (2 / 64)
        if boundary != -1:
            flux[:4] = flux[-4:] = 0.0
        expected = np.array([rho, momentum, momentum * temperature, flux])
        for a in coefficients:
            padded = np.pad(expected, ((0, 0), (4, 4)), mode="wrap" if boundary == -1 else "edge")
            expected = a[0] * padded[:, 4:68]
            for offset in range(1, 5):
                expected += a[offset] * (padded[:, 4 - offset : 68 - offset] + padded[:, 4 + offset : 68 + offset])
        actual = np.fromfile(tmp_path / f"restart_data/lustre_lso_stat_{step}.dat").reshape((-1, 64))
        np.testing.assert_allclose(actual[[1, 4, 7, 9]], expected, rtol=1e-11, atol=1e-11, err_msg=f"step={step}")
