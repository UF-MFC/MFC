# EL Particle Hemi-Sphere Shell Packing

This example initializes Euler-Lagrange solid particles with the in-code packing path instead of reading a particle `.dat` file.

The case uses:

- `particles_lagrange = "T"` to enable the EL solid particle solver.
- `lag_params%packing_flag = 2` for hemi-sphere shell packing.
- `lag_params%packing_size_distribution = 0` for constant-diameter particles.
- `lag_params%packing_volume_fraction` to estimate the generated particle count.
- `lag_params%packing_seed` for deterministic random placement.
- `lag_params%packing_shell_inner_radius` and `lag_params%packing_shell_outer_radius` for the shell geometry.

The generated particle centers are constrained so particle surfaces remain between the inner and outer shell radii and above the hemi-sphere plane. Particle-particle overlap is rejected through the same particle-cloud packing bridge used by the IBM packing path.

Run:

```bash
./mfc.sh validate examples/3D_lagrange_particle_hemi_shell_packing/case.py
./mfc.sh run examples/3D_lagrange_particle_hemi_shell_packing/case.py --clean --no-debug -- --vf 0.01 --seed 1 --steps 1
```

Useful case options:

```bash
--vf 0.01     # target packing volume fraction
--seed 1      # deterministic packing seed
--steps 1     # simulation stop/save step for quick checks
```
