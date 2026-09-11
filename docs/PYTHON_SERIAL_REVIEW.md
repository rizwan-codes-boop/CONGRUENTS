# Python serial preparation / native galaxy loops

Revision: 2026-09-11. Development branch:
`refactor/python-serial-galaxy-openmp`, based on the existing Week-2 work.

This is the preparation-stage review. The subsequent optional solver boundary
and current remaining work are documented in [WEEK3_REVIEW.md](WEEK3_REVIEW.md).

## Reviewed boundary

| Responsibility | Implementation now |
|---|---|
| Catalogue, Astropy input conversion, grid/temperature bounds | Python |
| Standalone ionisation diagnostic | Python |
| Radiation fields, dilution, integrated diagnostic densities | Python |
| IC emission/Gamma, bremsstrahlung and synchrotron table generation | Python, serial |
| Table interpolation, IC combination, NumPy storage and cache I/O | Python, serial |
| Galaxy properties previously computed in native galaxy loops | C/OpenMP over galaxies |
| Context allocation, status reporting and ABI declarations | Minimal C/ctypes plumbing |

Native helpers in `csrc/preparation.c` are reachable only from the galaxy-property
loop: gas surface density, velocity dispersion, stellar dispersion and dust
temperature, plus the embedded height/density/magnetic-field formulas.
These mirror the property work in `spectra.c` and `create_interp_objects.c`.
They must execute natively as part of the worker's loop body; they are not a
general-purpose C serial backend. There is exactly one OpenMP region in the
interface library. Its iteration index identifies a galaxy, not an energy bin.

Python prepares contiguous double rows [n,4] and output [n,10] before calling
`cg_galaxies`. The C thread-count variable is set before the pragma. Workers
read native inputs, use thread-local scalar temporaries and write disjoint
output rows. Nothing calls Python or converts units inside the parallel region.
Buffers live throughout the synchronous call; C retains no borrowed arrays.
The public header specifies column order and native units.

The SLUG methodology is retained: a shared library loaded with ctypes, explicit
argtypes/restype declarations, and straightforward Python wrappers. No Cython,
pybind11, Python callbacks, energy-bin OpenMP or parallel table generation.

## Numerical behaviour and compatibility

- Reference constants are copied into `constants.py` with explicit native
  values; Astropy remains for input conversion, not silent constant replacement.
- Python GK15 quadrature follows the legacy cubature rule, error estimator and
  scalar batch-refinement policy, at the original tolerances and evaluation
  limit. It is an attributed GPL adaptation; no C integration callback is used.
- SciPy supplies the modified Bessel function for synchrotron. NumPy supplies
  linear/bilinear interpolation. These compiled dependencies are ordinary
  Python scientific libraries, not a retained custom C table backend.
- CONGRUENTS-i dilution, emission/Gamma selection, diagnostic conventions and
  reversed temperature interpolation weights are preserved, including the
  explicit rejection of an unavailable upper temperature bracket.
- No calibration, output rescaling or scientific fixes were introduced.
- Standalone `Context.ionisation_loss` keeps its Python signature but is now
  serial. Thread count applies only to native galaxy processing.
- Library/package version is 0.2.0, ABI **2**. Removed C symbols include
  ionisation, table create/import/borrow/evaluate and IC combination.
  Old native clients must update; old shared libraries are rejected.
- NumPy table arrays are owned and read-only. Closing an object does not
  invalidate an already-retained NumPy view.
- Cache schema 2 fingerprints Python implementation files, NumPy/SciPy versions
  and the native library. Week-2 caches are not silently reused. Old files are
  left untouched.

## Validation

Run `make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python` in the astro
environment. All 11 catalogue rows are covered, not selected galaxies.
Local result: **23 tests passed**, including optional Astropy unit checks.

Independent legacy C reference comparisons use 32x32 tables for every IC
field/temperature and for both emission/Gamma families, plus BS and SY.
There are 26 planes for this catalogue. Combined IC tables are checked for
both kinds for all 11 galaxies. Reference sources do not call the new adapter.

Measured maximum relative differences on this Mac:

| Comparison | Maximum relative difference |
|---|---:|
| IC emission, 12,288 cells | 3.09e-12 |
| IC Gamma, 12,288 cells | 1.27e-11 |
| BS, 1,024 cells | 4.31e-15 |
| SY, 32 cells | 6.32e-14 |
| Combined IC emission, all 11 galaxies | 2.76e-13 |
| Combined IC Gamma, all 11 galaxies | 2.04e-12 |

The table regression tolerance is 5e-8 relative (absolute floor 1e-280), allowing
for the original 1e-8 integration criterion and floating-point differences,
not an observational-fit tolerance. Axes use 1e-12 relative tolerance.
Saved six/seven-digit diagnostics use 6e-7; the rounded production IC spot
fixture uses 1e-5. Ionisation uses 1e-14 against independent C evaluation.

Other checks cover analytical quadrature, all galaxy properties and radiation
diagnostics, one/four-thread equality, optional Astropy units, input validation,
native errors, cache corruption and reuse, table ownership and interpolation,
and the removal of serial-physics C entry points.

## Remaining acceptance gates

This updates the **already implemented preparation interface**, not the entire
CONGRUENTS solver. The electron solve, full galaxy-loop emission calculation,
final output handling and end-to-end spectra comparison remain future work.
Their native helper boundary must be reviewed against the galaxy-loop call graph.

Full 1000x1000 production Python precomputation has not been run or benchmarked.
Serial Python generation may be substantially slower; caches avoid repeating
it. The 32x32 comparisons and four-cell production fixture are not full-grid
or final-spectrum validation.

Local validation is macOS/arm64. CI is configured for a portable subset on
macOS and Linux, but its remote results and HPC execution remain unverified.
HPC acceptance requires the actual site's compiler, libraries and job allocation.
