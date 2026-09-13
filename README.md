# CONGRUENTS

A Python–C interface for modelling cosmic-ray transport and multiwavelength
emission from galaxies. The initial scientific target is reproduction of the
**improvised legacy code**, with its parameters, catalogue order and numerical
conventions preserved. The interface does not introduce a new calibration.

## Architecture

| Responsibility | Implementation |
|---|---|
| Catalogue, units, configuration, caching and output | Python |
| Radiation fields and IC, bremsstrahlung and synchrotron lookup tables | Python, serial |
| Normalisation, calorimetry, diffusion, primary/proton spectra and free-free inputs | Python |
| Galaxy properties, secondary injection, two-zone electron solves and source emission | C/OpenMP over galaxies |
| Internal absorption, EBL, observer conversion and diagnostic integration | Python |

Python loads shared libraries through `ctypes`, with explicit argument and return
types. It prepares contiguous float64 arrays and native descriptors before
entering OpenMP. Workers never call back into Python. There is no parallelism
over energy bins or lookup-table generation.

C contains executable loop bodies and their native helpers, not just OpenMP
pragmas. These include adaptive integration kernels, matrix assembly, the
energy-conserving electron solver, and IC/BS/SY/pion/neutrino emissivities.
Halo injection depends on the solved disc population. Removing further native
helpers requires a separately reviewed change to the staging or numerical
method; this is not claimed to be the smallest possible Python–C boundary.

Borrowed inputs remain Python-owned for the synchronous call and are neither
modified nor retained by C. Workers own their temporary matrices, splines and
accelerators and write disjoint galaxy rows. Calls to the solver backend are
serialized; galaxy work within each call is parallel. Preparation currently
covers the whole catalogue; memory-bounded batching is not implemented.

## Installation and build

Requirements: Python >=3.9, NumPy, SciPy, GNU GCC with OpenMP, GSL and cubature.
Astropy provides optional input-unit conversion. Reference constants remain
fixed; installing Astropy does not silently replace them. Matplotlib is needed
only for plotting. The GPL-2.0 licence and source attribution are retained.

```sh
python -m pip install -e '.[units]'
python -m pip install matplotlib  # optional plotting
```

The native code uses GNU nested functions, so Apple Clang cannot build the full
solver. The preparation library is version 0.2.0 / ABI 2; the solver uses ABI 3.
Rebuild shared libraries after interface changes.

### Existing local dependency installation

Run from this repository in the configured Python environment:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c shared solver PYTHON=python
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
```

The sibling directory supplies the compiler and libraries only. It is not an
input catalogue or a runtime scientific-output dependency.

### Linux and other installations

The dependency-free preparation subset can be built with:

```sh
make shared test-portable CC=gcc PYTHON=python
```

For the full solver, provide GSL and cubature installation prefixes:

```sh
make test CC=gcc PYTHON=python \
  CPPFLAGS="-I. -ICR_spectra -I/GSL_PREFIX/include -I/GSL_PREFIX/include/gsl -I/CUBATURE_PREFIX/include -DCONGRUENTS_USE_SYSTEM_OPENMP" \
  LDLIBS="/GSL_PREFIX/lib/libgsl.a /GSL_PREFIX/lib/libgslcblas.a -L/CUBATURE_PREFIX/lib -lcubature -lm"
```

Replace the placeholder paths. Static GSL archives must be built with
position-independent code (`--with-pic --disable-shared`) for use in the shared
solver library. The Linux recipe in
[the CI workflow](.github/workflows/interface.yml) builds a checksum-pinned GSL
2.8 from source and tests this configuration. Private static GSL linkage isolates
the solver's error handler from other libraries.

On macOS, use a GNU compiler such as Homebrew GCC 14 by supplying
`CC="$(brew --prefix gcc@14)/bin/gcc-14"`, together with compatible dependency
paths. On HPC, load the site's GCC/Python/library modules and limit
`--threads` to allocated CPUs. GNU nested callbacks may require executable
trampolines; target-site restrictions must be checked. HPC execution has not
yet been validated.

## Run the model

The bundled `input/cat_nt.txt` contains 11 model cases, including adjusted cases,
not 11 newly generated galaxies. Columns are redshift, stellar mass in solar
masses, effective radius in kpc, and star-formation rate in solar masses/year.
Preserve row order: the special index-10 proton-calorimetry factor of 0.1 is
part of the reference behaviour; primary-electron injection is not reduced.

```sh
# Fast regression example, not a converged physical prediction:
PYTHONPATH=src python -m congruents input/cat_nt.txt output/smoke \
  --threads 4 --tables 16 --cosmic-rays 16 --photons 8 --cells 16

# Production resolution: 1000x1000 tables, 1000 CR energies,
# 500 photon energies and 500 solver cells:
PYTHONPATH=src python -m congruents input/cat_nt.txt output/production \
  --threads 4 --cache data/python-cache

# Source-only Python API example:
PYTHONPATH=src python examples/source_spectra.py
```

Use a new output directory for each run; existing directories are rejected
before computation. Cache fingerprints include implementation/dependency
versions and the native library. Changed fingerprints regenerate tables rather
than silently reusing incompatible data. CANDELS integration is not implemented.

The complete `run(preparation, grid, threads=4)` API includes diagnostics.
The source-only `solve(...)` API leaves diagnostics disabled by default.
Results expose named read-only NumPy arrays; NPZ export does not use pickle.

## Outputs and units

Export preserves the numeric layout of 39 standard files from the improvised
legacy code, plus source NPZ, photon totals, optical depths and metadata.
Metadata records catalogue order, grids, precision, threads, units, fingerprints,
and the EBL input path/checksum.

| Output | Meaning and units |
|---|---|
| Energy axes | GeV; electron energy is total, CR/proton grid is kinetic |
| Source populations / injection | GeV^-1 / GeV^-1 s^-1 |
| Source emissivities | GeV^-1 s^-1 |
| Observer component spectra | E²-weighted flux, GeV cm^-2 s^-1 |
| `gal_data.txt` | h (pc), nH (cm^-3), B (G), gas dispersion (km/s), area (pc²), gas surface density (Msun/pc²), SFR surface density (Msun/yr/pc²), stellar surface density (Msun/pc²), dust T (K) |
| `Urad_Ub.txt` | Total radiation, magnetic, CMB, FIR, 3000 K, 4000 K, 7500 K, UV densities; eV/cm³ |
| `tau_loss/` | Loss times in seconds, galaxy rows and energy columns |
| `E_loss_nucrit.txt` | BS/SY/IC/ionisation/diffusion in disc, then BS/SY/IC/plasma/diffusion in halo; GeV/s at 1.49-GHz critical energies |
| `E_loss_leptons.txt` | Primary/secondary injected power (GeV/s), then primary and secondary SY/IC/BS disc, SY/IC/BS halo and escaped-injection fractions |
| `CR_specs.txt` | Five rows per galaxy: proton, primary/secondary disc, primary/secondary halo; GeV^-1 |
| `CR_specs_inj.txt` | Four rows per galaxy: primary/secondary disc and escaped halo injection; GeV^-1 s^-1 |
| `L_radio.txt` | Primary/secondary disc SY, primary/secondary halo SY, free-free; W/Hz at 1.49 GHz |
| `L_gamma.txt` | 0.1–100 GeV luminosity, GeV/s |

Zero tabulated loss rates produce positive infinite loss times. Other invalid
diagnostics raise errors. Native failures discard invalid galaxy results and
report affected indices rather than substituting accepted spectra.

`plot_figure9.py` reads files from `output/` relative to the working directory,
writes PNG/PDF there, and displays the figure. To plot `output/production`
without moving the data, use a temporary working directory whose `output`
symlink points to that directory, then run the script by its absolute path.
The plot uses six catalogue cases and the published panel arrangement; plotting
compatibility is not an independent test against digitised paper curves.

## Scientific provenance and preserved conventions

The improvised legacy code reconstructs missing constants, cosmology, math and
GSL-wrapper headers from call sites and units. Physical constants use CODATA
2018 / IAU conversion values. Flat GSL forwarding headers and the OpenMP
compatibility header remain necessary build infrastructure. The Python GK15
quadrature is an attributed GPL adaptation of the cubature rule and retains
its error/refinement policy and integration tolerances.

The reference already includes the following changes from the initially
recovered source; these were not introduced by the Python–C conversion:

- Radiation dilution denominator `2*pi` (Roth et al. 2023, Section 2.3).
- Free-free Gaunt-factor branch at 1 GHz and removal of the extra `1/f_e`
  after cancellation (equations 44, 45 and 51).
- Adjusted LMC proton-calorimetry scaling without tenfold primary-electron
  suppression (Section 5.2.2).
- Floating-point lookup validation using `fabs`, a tighter field-bound check,
  and CMB/FIR kernel selection by requested table type rather than field index.
- Per-galaxy interpolation accelerators and domain guards for pion and
  redshifted-spectrum interpolation.

Important reproduction conventions remain, even where they warrant scientific
review:

- Actual cosmology: H0=70 km/s/Mpc, Omega_m=0.3, Omega_Lambda=0.7, flat and
  without radiation. An older Planck-2018 source comment is not the calculation.
- Active EBL: Franceschini, despite an older Dominguez comment.
  `--ebl` explicitly selects an alternative. Energy is E/(1+z) in eV,
  capped at 1e15; bilinear extrapolation precedes nonnegative clipping.
- Internal absorption follows `spectra_funcs.h:tau_gg_gal_BW`, including its
  target-energy interval and scale-height path length.
- Observer conversion linearly resamples on E/(1+z), returns zero outside
  that domain and applies distance, E² and the original attenuation indexing.
  All 14 reference components, including neutrinos, receive attenuation.
  This is reproduced behaviour, not an endorsed neutrino prescription.
- Gamma luminosity includes internal absorption but excludes free-free,
  neutrinos and the full-calorimetry diagnostic. Photon totals include
  free-free and exclude neutrinos/full-calorimetry duplicates.
- Spectral integrals use adaptive log-energy quadrature of linearly
  interpolated spectra. Temperature interpolation retains the reference's
  reversed weights and rejects an unavailable upper bracket.
- Halo loss-time grids use ionisation; halo critical-energy diagnostics and
  the solver use plasma losses. Halo bremsstrahlung emission is absent.
- Default `legacy_table_precision=True` matches precompute-then-run:
  IC emission, BS and SY arrays undergo the six-decimal scientific-notation
  round trip. Gamma tables remain full precision because the reference loader
  rejects their saved bounds and regenerates them. This is not calibration.
  `False` selects an explicit full-precision diagnostic mode.
- Observer conversion requires z>0; independent distances for zero-redshift
  objects are not implemented.

Scientific corrections, changed constants, new catalogues and revised staging
must be proposed separately from interface maintenance.

## Validation

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
python tests/compare_outputs.py C_OUTPUT PYTHON_OUTPUT
make DEPENDENCY_ROOT=../CONGRUENTS-c test-runtime-sanitized
make DEPENDENCY_ROOT=../CONGRUENTS-c test-runtime-leaks  # macOS only
```

Tests independently compile the preserved C drivers with reduced resolutions
and extra exports; they do not use the new adapter to generate expected data.
Coverage includes all 11 cases, tables, transport, solved populations, emission,
observer spectra and diagnostics, units, invalid inputs, ownership, cache
reuse/corruption, error recovery, and one/four-thread equivalence.

Recorded validation on the merged implementation:

- 41 local macOS tests passed; latest pre-cleanup run: 52.538 s.
- Portable macOS/Ubuntu checks and the full Linux regression suite passed
  [CI](https://github.com/rizwan-codes-boop/CONGRUENTS/actions/runs/34705021553).
- Full macOS production resolution: all 39 files / 326,550 values passed against
  the improvised legacy code. Only two printed values differed, both in
  `CR_specs.txt`; maximum relative error 1.3698288624609998e-7.
- A cached repeat reproduced all exported numeric values. Initial run:
  92.67 min (25 generated tables, one reused); cached run: 32.60 min
  (26 hits). These are local observations, not controlled speedup measurements.
- Focused native runtime audit: 1,200 exercises passed, including sanitizer
  checks; the separate macOS leak audit reported zero leaks. This does not
  instrument the entire production solver or dependency internals.

Full exported comparisons require matching shapes, valid NaN/inf handling and
exact reference zeros, with rtol=3e-6 and atol=1e-280. Source comparisons use
2e-6; full small-grid tables use 5e-8. These are numerical tolerances, not
observational agreement percentages. No output rescaling was applied.

Regression fixtures are model outputs captured on 2026-09-06, not observations
or independently authenticated paper products. The `gal_data.txt` and
`Urad_Ub.txt` fixtures contain all 11 rows. Their SHA256 values are respectively
`11c524cb61bf99c46b8a3948503d8190e509d40c6d6c1799b6651b9c8826c5f6` and
`4576697395c1f8dd1521bc0808f5a9251dfebd0d9bfea6143fae243864ad1549`.
The IC JSON fixture records four cells at x indices 400/401 and y indices
200/201 from a 1000x1000 plane, with its source checksum and rounded coordinates.
Independent C comparisons supplement these historical fixtures.

Remaining validation: numerical convergence, full-resolution Linux/HPC
production, target-site portability and full-production peak-memory/leak
profiling. The earlier resource probe failed after successful export, so no
peak-memory measurement is claimed. Reproducing the improvised legacy code is
not independent validation of its physics or of published curves.

## Repository layout and maintenance

- `src/congruents/`: Python API, serial calculations and production workflow.
- `csrc/`: native galaxy workers, public ABI and checked solver runtime.
- `CR_spectra/`, root C sources/headers: improvised legacy code and build support;
  retained for native compilation and independent regression checks.
- `input/`: catalogue and three EBL datasets.
- `tests/`: automated checks, independent C drivers and numeric fixtures.
- `examples/source_spectra.py`: small source-spectrum API demonstration.
- `plot_figure9.py`: output-driven scientific figure.
- `LICENSE`: retained licence; preserve third-party attribution in source.

Generated binaries, caches, dependencies and outputs are ignored by Git.
Historical development notes and reports remain recoverable in Git history.

Standing maintenance requirements: keep the improvised legacy code as the
scientific reference; preserve its parameters and catalogue order; place work
outside OpenMP in Python and prepare native inputs before entry; parallelise
galaxies only; review any native-helper boundary change explicitly; validate all
cases and outputs; target macOS, Linux and HPC. Preserve unrelated user edits.
Pushes and merges require explicit user authorisation.
