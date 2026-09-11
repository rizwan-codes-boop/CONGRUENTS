# Week 3: galaxy solver and component emission

Branch: `feature/week3-solver-emission`. Reference: unchanged CONGRUENTS-i.

Week-3 implementation and local regression checks are complete. This is not
yet a validated full-production replacement; the acceptance gates below remain.

## Delivered

| Stage | Responsibility |
|---|---|
| Python setup | Catalogue, units, energy axes, serial lookup-table generation, IC combination, contiguous input/output arrays |
| Native transport galaxy loop | Proton calorimetry/diffusion, disc/halo electron diffusion, primary/secondary injection, steady-state protons |
| Native solver galaxy loop | Disc primary/secondary electron solve, escaped injection into halo, halo solve, component emission |
| Python result handling | Named read-only NumPy arrays, validation/errors, optional compressed NPZ export |

The native transport mirrors `spectra.c:242` onward and the injection work
inside `spectra.c:563` onward. The two-zone solver and emission follow
`spectra.c:629`, `spectra.c:652` and the following photon-energy loop.
There is no OpenMP over energy bins. The inner energy loops remain sequential
inside each native galaxy worker.

The baseline's catalogue-index-10 proton-calorimetry factor of 0.1 is retained.
Primary electron injection is not reduced. Catalogue order is therefore
scientifically significant; this API does not reinterpret a reordered catalogue.
No CANDELS input, output calibration or new physical parameters were added.

## Python/C boundary review

The existing lightweight library and ABI 2 remain unchanged. The optional
`libcongruents_solver` library has its own ABI 1 and requires GSL/cubature.
Both follow the SLUG ctypes loading pattern with explicit signatures.

Python prepares float64 contiguous arrays and native table descriptors before
calling C. References remain alive for the entire synchronous call; C retains
none after return. Each worker writes only its galaxy's output slices.
The native team size and ownership registries are prepared before OpenMP.
Galaxy-dependent solver matrices, splines and accelerator caches are created
inside native helpers, never via a callback to Python. This includes the helper
functions required to execute the loop body, not only the pragma itself.

The native call graph includes normalization/injection functions, radiative and
diffusive transition/loss kernels, the energy-conserving steady-state solver,
and IC/BS/SY/free-free/pion/neutrino emissivities. These functions execute only
as part of a galaxy worker. Serial lookup-table generators have no exported
entry point and are removed from the solver library by dead-code stripping.
Standalone preparation, table generation and combination remain Python.

`csrc/steady_state_native.h` is a reviewed copy of the energy-conserving part
of the original `CRe_steadystate.h`. Its nested LU helper is replaced with a
checked equivalent; allocation ownership is supplied by `solver_runtime.h`.
A test verifies the remaining equations against the original header verbatim.
The alternative number-conserving solver is not used or exported.

## API and units

`solve(preparation, SolverGrid(...), threads=4)` returns `SolverResult`.
Default grids retain production counts: 1000 cosmic-ray energies, 500 photon
energies and 500 solver cells; table resolution is controlled by Preparation.
The example and tests explicitly use smaller grids and are not converged spectra.

| Result | Shape / units |
|---|---|
| Kinetic/electron/photon axes | GeV; electron energy is total, not kinetic |
| fcal | [galaxy, CR energy], dimensionless |
| Proton, disc-electron, halo-electron diffusion | [galaxy, CR energy], cm²/s |
| Primary/secondary disc injection | [galaxy, CR energy], GeV⁻¹ s⁻¹ |
| Proton steady state | [galaxy, CR energy], GeV⁻¹ |
| Four electron populations | [galaxy, CR energy], GeV⁻¹ |
| IC, BS, SY, free-free, pion and neutrino emission | [galaxy, photon energy], GeV⁻¹ s⁻¹ |
| Free-free optical depth | [galaxy, photon energy], dimensionless |

Disc IC/BS/SY include the reference free-free attenuation. Halo IC/SY do not.
The reference does not emit a halo-BS component in this stage. Pion emission
includes both the ordinary and full-calorimetry diagnostic spectra.
These are source components, not observer-frame fluxes or summed E² spectra.
`result.save(path)` writes named arrays in NPZ format without pickle. It does
not claim compatibility with the legacy output-file layout.
Metadata records the catalogue order, grids, precision convention, units,
thread count and preparation/native-library fingerprints.

## Table precision: a necessary compatibility detail

The default `legacy_table_precision=True` reproduces the normal
`run-precompute` then `run` workflow:

1. IC emission, BS and SY planes/axes are rounded through the equivalent of
   the legacy `%le` writer (six decimal places in scientific notation).
2. Gamma tables stay full precision. The reference loader checks their axes
   against emission limits, rejects the stored Gamma files and regenerates
   them. See `data_calc.h:200` and the corresponding 3-D check.
3. Temperature axes follow the corresponding stored/generated convention.

Rounding happens on independent temporary Python arrays; the preparation
cache and original files are untouched. It is not a fitted normalization.
Without matching this convention, the rounded table endpoints can exclude
the first/last photon-grid coordinate and give zeros where full precision does
not. Tests cover these endpoints, rather than discarding them.

`legacy_table_precision=False` is an explicit full-precision diagnostic mode.
It is tested against a separately generated C reference whose table writer
uses 17 digits. It is not the default production-compatibility mode.
A cold legacy `spectra` run with missing tables may itself retain freshly
computed, unrounded planes; the default deliberately targets precompute + run.

## Error handling and lifetimes

- Python validates grid sizes, threads, preparation lifetime and library ABI.
- C validates axes, dimensions, values and required buffers before workers run.
- Each worker tracks allocations and releases them on success or failure.
- Non-finite integrands return through cubature first so its workspace can be
  released. The adapter then unwinds the worker, not the Python process.
  A runtime audit found that cubature's early callback-error path leaked memory.
  The wrapper now records the error and supplies temporary zero callback values
  so cubature can finish its normal cleanup. The integration result is marked
  invalid and discarded; these temporary values never become accepted spectra.
- LU failures/non-finite results become per-galaxy errors. Invalid spline
  evaluation returns a numerical error instead of the reference's abort.
- The adapter discards outputs if any galaxy fails and reports its indices.
- GSL symbols are private to this backend. Its handler is disabled and restored
  under a batch guard; the legacy abort handler is never installed.
- Calls to the optional backend are serialized; OpenMP parallelises the galaxies
  within a call. This avoids conflicting backend-global handler changes.

An 8-cell trial gives non-finite halo integrands for several reference galaxies.
The test checks that this fails safely and a subsequent 16-cell run succeeds.
Neither a successful 16-cell run nor clipping a bad result constitutes a
convergence test. The adapter does not silently replace non-finite values.

## Reproduce the checks

In the astro environment:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test-week3 PYTHON=python
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
make DEPENDENCY_ROOT=../CONGRUENTS-c test-runtime-sanitized
make DEPENDENCY_ROOT=../CONGRUENTS-c test-runtime-leaks # macOS only
PYTHONPATH=src python examples/week3.py
```

The reference build uses the original production drivers/headers independently,
not the new C adapter or Python equations. The build script changes resolution,
adds raw-output statements and stops before observer-frame postprocessing.
The full-precision comparison additionally changes only the table writer.
Reference runs use temporary directories, never existing production outputs.

All 11 galaxies are compared for all 26 exposed arrays: seven transport/proton,
four electron and fifteen emission/optical-depth arrays. The 16-cell/16-CR-energy/
8-photon-energy comparison checks 3,256 entries per precision convention, using
16x16 tables. A second legacy-precision comparison uses 32 cells, 32 CR energies,
16 photon energies and 32x32 tables, checking another 6,512 entries. These are
regression grids, not evidence of production convergence.
The tolerance is 2e-6 relative with a 1e-280
absolute floor, allowing for six-digit reference output and numerical
integration; exact zero endpoints are included. One/four-thread results must
be bitwise equal. Other tests cover export, read-only outputs, raw invalid
arguments, ABI/symbol isolation, solver-equation fidelity and failure recovery.

The final local suite passed all 31 Python tests. A separate native runtime
harness passed 1,200 allocation, interpolation, LU, non-finite-integrand and
failure/recovery exercises, including injected allocation failures. The same
harness passed AddressSanitizer and UndefinedBehaviorSanitizer; macOS `leaks`
reported zero leaks and zero leaked bytes. macOS AddressSanitizer runs with
`detect_leaks=0`, so the separate leak audit is important. These focused checks
do not instrument the entire production solver or third-party library internals.

## Remaining Week-4 / production acceptance

- Full 1000x1000 Python table generation and 500-cell/all-galaxy production
  comparison, convergence and performance measurements have not been run.
- EBL/gamma-gamma attenuation, observer distances/redshifting, total spectra,
  the complete loss-time/energy-budget/radio diagnostics and legacy output
  layout still need Python orchestration and all-output validation.
- Fresh Linux/HPC execution and full-production solver sanitizer coverage remain outstanding.
  Local validation is macOS/arm64. Existing CI only covers the portable
  preparation subset, not the new GSL-dependent solver.
- No new merge or push is authorized by the Week-3 implementation request.

For Linux/HPC, supply the site's GCC, GSL and cubature include/library flags,
including the flat GSL include directory required by the legacy headers:

```sh
make solver CC=gcc \
  CPPFLAGS="-I. -ICR_spectra -I/GSL_PREFIX/include/gsl -I/CUBATURE_PREFIX/include" \
  LDLIBS="-L/GSL_PREFIX/lib -L/CUBATURE_PREFIX/lib -lgsl -lgslcblas -lcubature -lm"
```

Replace the prefix placeholders and load the required modules/runtime paths.
Use no more threads than the job allocation. GNU nested callbacks are retained
from the reference; site executable-trampoline restrictions must be checked.
