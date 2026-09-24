@page lso-filter-testing LSO Filter Verification Plan

# LSO Filter Verification Plan

This page defines the acceptance tests for the LSO filtered statistical products and
Euler-Lagrange closure output. A successful build or a visually smooth field is not a
correctness test. Production use requires the numerical, MPI, GPU, and physics checks below.

## Supported scope

The simulation forms each nonlinear statistical product from the unfiltered state and then
filters the product. Post-process reads those simulation products. An optional post-process
pass can widen the filter by filtering the products and the gas mask again.

The 11 product blocks currently require one fluid. They use `lso_R_gas` for ideal and stiffened gases and the common EOS
temperature interface for state-dependent families. Closure reconstruction is currently
limited to one calorically perfect ideal or stiffened gas. Chemistry, state-dependent
EOS closures, and stretched grids remain outside the closure acceptance scope. Serial LSO
output is implemented. Simulation downsampling uses coarse-grid MPI views and requires
divisibility in every active grid direction. Reduced-grid post-processing requires
`parallel_io=T`, `file_per_process=F`, at least two coarse cells in each active direction,
and a valid MPI decomposition of the coarse grid. It rejects legacy `down_sample=T`
and spatial boundary-condition files. Other layouts remain outside this reader's scope.
The simulation's coarse-cell-center interpolation is unchanged.

## Audit against upstream, 2026-09-20

Reference: MFlowCode/MFC commit `8551f17ab4866674ac88ffc4032dd83a0050e7c5`.
The fork worktree includes an uncommitted EOS refactor integration. The central
conversion, HLLC, Python EOS and reactive-burn files match that upstream snapshot;
this is provenance, not experimental validation. The initial audit left upstream
goldens unchanged; the later user-authorized IBM candidate is recorded below.
The two added LSO goldens cover filtered conserved output, not all
statistical products or closure quantities.

The post-process closure loop reused the final cell's mask weight. A 32 by 32,
two-step, periodic viscous IBM test independently reconstructs both components of
W_tau_u from the saved statistical fields and each cell's saved gas fraction.
It also checks that the particle fraction and gas mask sum to one.
The fix reloads the weight in the closure loop and
does not change the flow equations, time integration or filter coefficients.
The stiffened-gas case also checks Favre temperature against the saved energy and
kinetic-energy moment. Its thermal offset is the physical stiffening pressure,
obtained with the central EOS coefficient helper, not the stored energy coefficient.

Run this executable regression explicitly after building all three targets:

```shell
MFC_TEST_LSO_OUTPUT=1 PYTEST_ADDOPTS='-k test_viscous_closure_uses_each_cells_mask' ./mfc.sh lint
```

Each case has a 180-second execution timeout, uses a fresh temporary directory, checks all
21 statistical components and 17 closure components for finiteness, and checks the
two viscous-work components and temperature against their algebraic definitions for
ideal and stiffened gas with periodic and extrapolated boundaries, at stride 1 and 2.
It does not validate all remaining closure formulas, 3D IBM or GPU execution.
Ordinary lint skips it because lint must not require compiled solver executables.

The nonperiodic variants initially failed the phase-partition check in 636 of 1024
cells, with maximum error 0.16853033. Statistical convolution froze a boundary
band while flow and mask convolution updated every interior cell. Removing that
band restriction restores the same linear operator and ghost extrapolation for
all three. All eight variants pass, including phase partition at `1e-10`
absolute tolerance. No filter coefficients, flow equations or downsampling changed.

### Reduced-grid reader and GPU scoping correction

The old post-process reader used fine-grid dimensions for coarse filtered files.
A stride-2 run produced twice the expected extent in each active direction; the
density output contained adjacent momentum-field data despite a successful exit.
Post-process now decomposes the coarse grid, reads matching field offsets, and
uses every stride-th boundary from the original grid file. The initial unfiltered
snapshot is interpolated at coarse-cell centers using sparse MPI reads and
coarse-sized temporary buffers. Categorical IBM markers use the nearest fine cell,
with the lower index chosen for an even-stride tie; they are not averaged.

An independent visualization defect appeared when reading the initial and filtered
Silo snapshots consecutively: the field-layout cache reused initial field locations
after LSO fields were inserted. It now caches by rank and file-object set, preserving
reuse within a stable layout. The executable reader test covers this transition.

The 1D statistical stress used `2*mu*du/dx`, while the resolved closure used
`4*mu*du/dx/3`. Both statistical implementations now use the latter, consistently
with the multidimensional deviatoric-stress convention. A 64-cell linear-velocity
case (`u=0.1*x`, `mu=0.01`) checks interior stress against `0.0013333333333333333`
and viscous subfilter stress against zero. The old code gave approximately
`0.002` and `0.0006666667`, respectively. This changes diagnostic statistics,
not the solver's viscous fluxes.

The algebraic and 1D/2D/3D gradient kernels explicitly privatize per-cell scratch
variables with `GPU_PARALLEL_LOOP`. No filter arithmetic, coefficients, time-step
control or IBM state correction changes with this GPU patch. Macro expansion was
inspected; accelerator compilation, device execution and performance are still
required on supported hardware, not inferred from a Mac CPU pass.

Run the bounded executable suite after building all three targets:

```shell
MFC_TEST_LSO_OUTPUT=1 PYTEST_ADDOPTS='-k test_lso_closure_output' ./mfc.sh lint
```

It contains ten IBM closure/mask/marker checks, six reader cases (1D stride 2
on one rank, 1D stride 3 on two ranks, and 2D/3D stride 2 on two ranks), and the
linear-velocity stress check, and two saved-state temperature-product checks.
The additional reader and IBM cases cover adaptive and constant CFL without
fixed-step time parameters. Reader checks compare every saved conservative
field with raw data at initial and filtered times, including MPI interface cells,
and verify coarse coordinates. The multidimensional fixtures vary in each active
direction. Each subprocess has a 180-second timeout; budget 10 minutes for this
laptop test batch with four build jobs. Do not regenerate goldens for these fixes.

The 3D reader case runs with bounds checking. Its first debug run found that
WENO3 allocated three halo cells while the LSO stencil accessed four. Simulation
now raises the halo allocation and active-direction bounds to at least four cells
when LSO is enabled, before MPI buffers and device data are initialized. Non-LSO
allocation is unchanged. The debug case is mandatory in this opt-in suite because
release-only raw-file comparisons did not expose the out-of-bounds access.

Final local checks on 2026-09-20: all 13 executable regressions passed in 41.17 s,
including the two-rank 3D debug case. The eight UUIDs listed below passed again in
17.3 s without regenerating goldens. All three CPU targets built successfully.
These are targeted double-precision checks, not the complete 777-case suite or
validation of every EOS, closure, particle model, compiler and GPU backend.

Files changed for this reader/stress/GPU correction:

- Simulation: `m_global_parameters.fpp`, `m_lso_filter.fpp`.
- Post-process: `m_data_input.f90`, `m_start_up.fpp`, `m_lso_pp_filter.fpp`.
- Toolchain: `case_validator.py`, `test_case_validator.py`,
  `test_lso_closure_output.py`, `viz/silo_reader.py`.
- Documentation: `case.md` and this page.

The sparse MPI-view helper is shared by initial-field and IBM-marker reads. No
new dependency or configuration option is added; no further safe deletion was
identified within this patch. Earlier uncommitted integration changes are separate.

### Saved-state diagnostics and CFL indexing

The next audit reproduced a stale RK-stage pressure in the temperature products:
the 64-cell shock case differed from a saved-energy reconstruction by 0.0094781
at step 4. Statistics now reconstruct primitives with the central conversion in
the existing filter workspace, before reusing that workspace for filtering.
One additional temperature field holds the saved-state values and exchanges
their MPI halos for heat-flux gradients. Live RK primitives and temperatures,
Dan's IBM corrections, EOS formulas and filter/coarsening coefficients are unchanged.
This work occurs only when statistical output is saved.

The temperature regression independently forms rho, momentum, rho*u*T and
Fourier heat flux from raw conserved snapshots, then applies the recorded FIR
coefficients. It covers one-rank extrapolated and two-rank periodic/conducting
cases at steps 1 and 4, with a 1e-11 comparison tolerance.

CFL post-processing now uses n_start and unit save-index spacing. Shared IBM
marker writes use the save index directly in CFL mode. Both CFL modes previously
requested a missing initial filtered file and left the expected marker offsets
empty. The new tests exercise both failures; fixed-step paths remain covered.

Filtered conserved fields and the IBM gas mask are copied from device to host
before CPU sampling or writing. Previously that transfer copied the original
conserved fields instead. Separate-memory accelerator validation remains required.

Local verification on 2026-09-20 (M4 Pro, gfortran 15.2, CPU double): all three
targets built, precheck passed 7/7, and the 19 executable checks passed in 53.61 s.
The eight selected goldens passed in 26.1 s without regeneration. Against the
pre-fix audit run, all five raw and four filtered conserved snapshots are
byte-identical; the step-4 rho*u*T error fell from 0.0094781 to 1.07e-14.
These checks do not certify chemistry, mixed precision or accelerator execution.

This correction changes simulation's `m_lso_filter.fpp`, `m_start_up.fpp` and
`m_data_output.fpp`; post-process's `m_start_up.fpp`;
`toolchain/mfc/test_lso_closure_output.py`; and this page. It adds no settings
or dependencies. No further safe deletion was identified within this correction.

Remaining release blockers:

- Multi-fluid statistics used first-fluid partial density with mixture momentum;
  validation now rejects that unsupported combination. Filtered flow output alone
  remains available for multiple fluids.
- Solid Euler-Lagrange particles and Lagrangian bubbles are different paths. The
  local MG bubble golden does not validate solid-particle drag, heat transfer,
  collision conservation or restart migration.
- In `s_compute_lso_stat_fields`, particle fraction and particle velocity come
  from IBM markers only. With solid EL particles alone, that branch sets both
  products to zero. Validation now rejects `particles_lagrange=T` with `lso_stat_wrt=T`
  rather than writing an incomplete solid-EL particle closure dataset.
- Moving-IBM same-temperature density reconstruction has four Newton iterations
  without a final residual check. Its supported state range needs independent tests.
- Moving-IBM divergence diagnostics perform global gathers each step; assess their
  production cost and resolve cross-rank consistency before removing diagnostics.
- All 11 products need manufactured-field checks before and after widening, plus
  refinement and MPI/backend comparisons. Smooth contours and smoke passes do not
  establish physical correctness.
- The coarse second-stage cascade now also filters the decimated statistical products
  with the same operator as conserved variables and the mask. This restores matching
  widths without changing stride sampling; broader manufactured-product tests remain necessary.
- The viscous-product construction still zeros stress, heat flux and viscous-power
  products in a physical-boundary band before convolution. The phase-partition fix
  does not validate those gradients. Nonuniform boundary-gradient tests and a
  reviewed boundary treatment remain necessary.

Verification and physical validation are separate gates; see the
[NASA CFD verification/validation tutorial](https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html)
and its [grid-convergence guidance](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).
Use independent analytic inert-flow cases and measured particle-flow data with
uncertainty estimates, not newly generated fork output as its own reference.

### Initial local checks, before the IBM candidate update

macOS CPU, MPI-enabled, double precision, release build, no case optimization:

| Check | Result |
| --- | --- |
| Format and all seven precheck gates | Passed |
| Build pre_process, simulation, post_process | Passed |
| New periodic ideal/stiffened-gas closure regressions | 2 passed in 5.28 s; both failed against the old executable |
| Boundary-filter follow-up, periodic and extrapolated | 4 passed in 8.74 s; both extrapolated variants failed before this fix |
| Selected golden batch, including post-process | 7 passed, 1 failed in 26.4 s |
| Existing upstream golden content | Unchanged; no regeneration performed |

The passing UUIDs were `5AC2F65D`, `ED75D01D`, `471270BB`, `120043F6`,
`9521F2F8`, `FCBFF479`, and `AA49A8BC`. Moving-cylinder case `127A967A`
failed against the original reference: maximum failing absolute error 0.05641878 in conserved energy;
maximum failing relative error about 0.2171 in momentum (absolute error 0.00464023).
The comparison tolerance was 0.001 absolute and relative. Do not refresh its golden
to hide this discrepancy. This is a selected regression batch, not the full suite.

```shell
./mfc.sh test -j 4 -a --no-build --only 5AC2F65D ED75D01D 471270BB 120043F6 9521F2F8 FCBFF479 AA49A8BC 127A967A
```

Changed files in this audit pass (excluding the pre-existing staged integration):

- `src/post_process/m_lso_pp_filter.fpp`: per-cell normalization and thermal offset.
- `src/simulation/m_lso_filter.fpp`: matching statistical convolution extent at physical boundaries.
- `toolchain/mfc/case_validator.py`: single-fluid statistics guard and EL input documentation.
- `toolchain/mfc/test_case_validator.py`: unsupported mixture regression.
- `toolchain/mfc/test_lso_closure_output.py`: executable algebraic regression.
- `.github/scripts/monitor_slurm_job.sh`: portable file-size query.
- `toolchain/mfc/test_submit_requeue.py`: normalize the fixture's padded line count.
- `docs/documentation/case.md`: solver choices and LSO input/output descriptions.
- `docs/documentation/lso-filter-testing.md`: audit evidence and release blockers.

No solver abstractions or runtime options were added. Temporary IBM global-gather
diagnostics remain a candidate for later removal after cross-rank consistency is verified.

### Moving-cylinder discrepancy isolated

A fresh checkout of upstream `8551f17a`, built with GNU Fortran 15.2.0 on the
same Mac, passes `127A967A`; the rebuilt fork fails it. Both use CPU release,
MPI, double precision, no fast math and no case optimization. The diagnostic
comparison uses one rank, the same fixed `dt = 6e-6`, and identical generated
simulation inputs, with every timestep saved and body-state output enabled.

Initial fields and body states agree. The first saved difference is at step 1:
maximum density difference `2.9713424075e-4`, energy difference `7.4430015545e-4`.
This case does not exercise adaptive CFL or a state-dependent EOS.

In an isolated upstream checkout, an eight-insertion, one-deletion diagnostic
patch reproduces the ideal-gas specialization of the fork's moving-wall state
changes: pressure-consistent density, matching conserved energy, and synchronized
primitive velocity. With that patch, all 255 conserved-field text files and all
51 binary body-state files are byte-identical to the fork over steps 0 through 50.
The diagnostic patch is not a general-EOS implementation and was not applied to
the fork. The density/energy changes trace to `df5fc461d`; primitive-velocity
synchronization traces to `800686c09`.

This establishes causality for this regression, not physical validation of the
revised boundary treatment. Retaining these changes and demanding bit-identical
agreement with the old upstream result are incompatible for this case. Preserve
the original upstream reference in its pinned commit. The user selected preservation of Dan's corrections
on 2026-09-20; independent validation remains required before declaring the
revised method production-ready. Do not treat the current mismatch as an
unexplained platform failure, report it as a passing regression, or refresh its
reference to hide the difference.

On 2026-09-20, both requested branches were fetched from `danieljvickers/MFC-Dan`:
`debug-ibm-stability` at `033e680c2de9b41f03b79ccd15918ed7e8a0e757` and
`ibm-stability-integration` at `c148b2cc9e6abf2420bbf1ae743ffd1aa8abfacd`.
Both contain the same `127A967A` golden as upstream and the fork before regeneration (Git blob
`461c051fce852bf880aea8d198f5b7f2ef9a6c04`). Its metadata records generation on
2025-12-04 from a dirty `forces-via-volume-integrals` tree at `cfe671402d`.
Neither branch supplies an updated reference for the preserved corrections.
Rebuilding all three targets and rerunning `./mfc.sh test -j 4 -a --only 127A967A`
again fails against that identical reference (4.4 seconds of testing; maximum
failing energy difference `0.05641877816683`). No golden was replaced during that comparison.

### User-authorized moving-IBM candidate golden

The user subsequently authorized a new fork reference. Only `127A967A` was
regenerated with `./mfc.sh test -j 1 -a --no-build --only 127A967A --generate`.
Its new Git blob is `0e58b1edba05cd0d13f95516066bc456fd38744e`.
The original is recoverable from the pinned upstream commit above; no other
existing golden was regenerated. The generated metadata accurately identifies
the current dirty worktree, compiler and configuration. This is a provisional
fork regression reference, not a golden obtained from Dan or a release validation.

An ordinary run without regeneration passed all eight UUIDs in the recorded
batch above in 25.9 seconds, including pre-process, simulation and post-process.
All ten initial/final conserved-field entries in the new golden exactly match
the earlier stepwise diagnostic. Across its 51 saved states, all conserved fields
are finite, density is positive (minimum 0.9526646882), and independently
reconstructed ideal-gas pressure is positive (minimum 0.9350078995).

The largest step-50 energy difference from the pristine upstream diagnostic lies
inside the body (distance from its center is 0.6152 radii). For the 672 cell centers
outside the body in both runs, the maximum absolute energy difference is
0.0084698929 and sum(abs(delta E))/sum(abs(E_upstream)) is 0.0000496210.
These checks explain why a local golden failure can coexist with a small bulk
field difference; they do not validate wall force, conservation or convergence.
Independent moving-wall benchmarks and MPI/backend checks remain release gates.

## Numerical invariants

Every implementation change must verify:

1. Constant preservation: filtering a constant field returns the same field to roundoff.
2. Kernel normalization: each pass satisfies
   `a0 + 2*(a1 + a2 + a3 + a4) = 1`.
3. Design convergence: filter generation fails when the requested transfer-function error
   is not reached within `lso_max_passes`.
4. Phase partition: `phi_p + filter(g) = 1` to the output precision.
5. Product ordering: `rho_uu` equals `filter(g*rho*u*u)`, not a product formed from filtered
   density and momentum.
6. Mask normalization: filtered conserved variables equal
   `filter(g*q)/max(filter(g), lso_w_floor)`.
7. Widening equivalence: the post-process pass gives
   `filter2(filter1(g*X))` for every product and uses
   `filter2(filter1(g))` for normalization.
8. Uniform-flow limit: `R_sg`, `Q_T`, `E_ku`, `W_tau_u`, `R_mu_sg`, and `R_lam_sg` are zero
   to roundoff away from physical boundaries.
9. Manufactured nonuniform field: each closure agrees with an independent NumPy evaluation
   using the saved filter coefficients.

The negative version of the product-ordering test is required: constructing `rho_uu` after
filtering must fail the manufactured-field comparison.

## Functional cases

Run these cases with bounds checking:

```shell
./mfc.sh build -j 4 --debug
./mfc.sh run examples/1D_lso_mach3_shocktube/case.py -n 2 --debug -j 4
./mfc.sh run examples/1D_lso_mach3_shocktube/case.py -n 2 -t post_process --debug -j 4
```

The one-dimensional case verifies nonperiodic shock-tube boundaries, nonzero sub-filter
stress near the shock, finite closure output, and post-process reading of `lso_stat` data.
The initial restart has no filtered companion file, so the initial database must contain
the ordinary flow fields and omit LSO statistics and closures. Every later saved step must
contain exactly 11 statistical blocks and the dimension-appropriate closure count.

The reduced three-dimensional IBM case must contain at least one resolved sphere and use
periodic transverse boundaries. It verifies `phi_p`, stationary `phi_p*u_p = 0`, mask
normalization, and all 3D tensor components. A moving-body variant must verify that the
saved in-situ particle velocity, rather than the case-file initial velocity, reaches
post-process output.

## MPI tests

Run the manufactured and reduced IBM cases with one, two, and four ranks. Compare every
filtered conserved field, mask, statistical product, and closure after assembling the
global arrays. The expected result is bitwise equality when the decomposition is unchanged
by I/O packing, otherwise agreement within the configured storage-precision tolerance.

Use a debug build for the two-rank shock-tube test. It specifically protects the global
QOI-layout initialization, compact MPI reads into buffered fields, four-cell filter halos,
and closure output split into communication-sized chunks.

Include decompositions that cut through the shock and through an immersed body. Verify
periodic wrapping in y and z and nonperiodic extrapolation in x independently.

## GPU and compiler matrix

The production gate consists of:

| Lane | Required result |
| --- | --- |
| gfortran CPU debug and release | Build and numerical checks pass |
| nvfortran OpenACC | Build and numerical checks pass |
| nvfortran OpenMP target | Build and numerical checks pass |
| Cray OpenACC and OpenMP target | Build and numerical checks pass |
| AMD flang OpenMP target | Build and numerical checks pass |
| Intel ifx CPU | Build and numerical checks pass |
| Case optimization, CPU and available GPU lanes | Same supported-case results |

CPU success is not evidence of GPU race freedom. Compare device output with the CPU
reference for every statistical field, including the viscous products.

## Validation tests

Unit tests must reject:

- statistics without filtered output;
- statistical output without parallel I/O;
- nonpositive `filter_sigma`;
- stretched grids;
- non-divisible LSO downsampling;
- closures without statistics;
- closures for chemistry, multiple fluids, state-dependent EOS families, or nonpositive
  `lso_R_gas`.

## Release gate

Run, in order:

```shell
./mfc.sh format -j 8
./mfc.sh precheck -j 8
./mfc.sh build -j 8 --debug
./mfc.sh test -j 8 --only <LSO test UUIDs>
```

No golden may be generated until the independent numerical checks pass. Generate only the
LSO goldens from a clean checkout of the reviewed commit. Record compiler, backend, rank
count, precision, case-optimization state, and the maximum error for each invariant.

<div style='text-align:center; font-size:0.75rem; color:#888; padding:16px 0 0;'>Page last updated: 2026-09-20</div>
