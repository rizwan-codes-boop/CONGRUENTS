# Week 4: observer workflow and production acceptance

Branch: `feature/week4-observer-workflow`. Week 3 was published and merged into
main at `73d64d5`. Week-4 implementation and local macOS production-reproduction
checks are complete. Cross-platform execution and convergence certification
remain unverified, as detailed below.

Implementation commit: `98ce3e5`. The existing `plot_figure9.py` also ran
unchanged against the full-size Python output, producing PNG/PDF figures.
This checks plotting compatibility; it is not independent paper validation.

## Delivered workflow

`python -m congruents CATALOGUE NEW_OUTPUT_DIRECTORY` runs the source solver,
attenuation, observer transformation and diagnostic export. It defaults to
1000x1000 tables, 1000 CR energies, 500 photon energies and 500 solver cells.
Use `--cache` for reusable, fingerprinted Python tables and `--threads` for the
galaxy OpenMP team. Existing output directories are rejected before computation.

| Component | Implementation / contract |
|---|---|
| Source preparation | Python, unchanged Week-3 precomputation boundary |
| Galaxy solves and source emission | Existing native workers; optional 1.49-GHz SY query added |
| EBL | Python parser, bilinear extrapolation, energy cap and clipping |
| Internal absorption | Python Breit–Wheeler kernel and adaptive log-energy quadrature |
| Observer conversion | Python distances, linear redshift resampling and attenuation |
| Diagnostics | Python loss-time integrals, energy budgets, CR arrays, radiation densities and radio conversion |
| Export | 39 standard legacy files, plus source NPZ, photon total, optical depths and provenance |

Solver ABI is now **3**; the preparation ABI remains 2. Rebuild with `make solver`.
The new optional native buffer contains four *raw* synchrotron emissivities at
1.49 GHz, evaluated from the solved disc/halo spectra. Python supplies free-free
emission/attenuation and unit conversion. No diagnostic quadrature, table
generation or serial output processing has been moved into C. All native
parallelism remains over galaxies, without Python callbacks.

The source `solve(...)` API leaves diagnostics disabled by default; the complete
`run(...)` API enables them. Diagnostic arrays are available as
`result.source.diagnostics`. Zero tabulated loss rates give positive infinite
loss times, not silently substituted finite values. Other invalid diagnostics
raise errors.

### Diagnostic columns and units

| File/group | Columns / units |
|---|---|
| gal_data | h (pc), nH (cm^-3), B (G), gas dispersion (km/s), area (pc²), gas surface density (Msun/pc²), SFR surface density (Msun/yr/pc²), stellar surface density (Msun/pc²), dust T (K) |
| Urad_Ub | total radiation, magnetic, CMB, FIR, 3000 K, 4000 K, 7500 K, UV; eV/cm³ |
| tau_loss/* | Loss times in seconds; rows are galaxies, columns are total electron energies (proton diffusion uses the kinetic CR grid) |
| E_loss_nucrit | BS/SY/IC/ionisation/diffusion in disc, then BS/SY/IC/plasma/diffusion in halo; GeV/s at each zone's 1.49-GHz critical energy |
| E_loss_leptons | Primary/secondary injected power (GeV/s), followed by primary then secondary SY/IC/BS disc, SY/IC/BS halo and escaped-injection fractions; halo BS is zero |
| CR_specs | Five rows per galaxy: proton, primary/secondary disc, primary/secondary halo; GeV^-1 |
| CR_specs_inj | Four rows per galaxy: primary/secondary disc, primary/secondary escaped halo injection; GeV^-1 s^-1 |
| L_radio | Primary/secondary disc SY, primary/secondary halo SY, free-free; W/Hz at 1.49 GHz |
| L_gamma | 0.1–100 GeV luminosity; GeV/s |

Numeric shapes, ordering and precision match the legacy files. Export header
wording is not required to be byte-identical; the comparisons check numeric data.

## Reference conventions deliberately preserved

- Cosmology uses the actual `cosmo_params.h` values: H0=70 km/s/Mpc, Om=0.3,
  Ode=0.7, flat and without radiation. The Planck-2018 comment is not the model.
- The active EBL input is **Franceschini**, despite a Dominguez comment in
  `spectra.c`. Users may explicitly supply another table with `--ebl`.
  Its absolute path and SHA256 are recorded in metadata.
- EBL energy is E/(1+z) in eV, capped at 1e15; interpolation extrapolates
  bilinearly outside the table before clipping optical depth to nonnegative.
- Internal optical depth follows `spectra_funcs.h:tau_gg_gal_BW`, including
  the original target-energy integration interval and scale-height path length.
- Observer conversion follows `spectra.c:1003`, `spectra.c:1028` and
  `file_io.h:write_2D_spec_file`: linear interpolation on E/(1+z), zero outside
  its contracted domain, multiply by both attenuation arrays, distance factor
  and E². Optical depths retain the writer's original array indexing.
- The legacy writer attenuates **all 14 components, including neutrinos**.
  This is retained for reproduction, not endorsed as a physical prescription.
- `L_gamma` integrates 0.1–100 GeV with internal absorption only; it excludes
  free-free, neutrinos and the full-calorimetry diagnostic. The separately
  exported photon total includes free-free and excludes neutrinos/fcal1.
- Spectral integration follows the original adaptive log-energy quadrature of
  a *linearly* interpolated spectrum, not log-log interpolation or a trapezoid.
- The halo loss-time grid uses ionisation, while its critical-energy diagnostic
  and solver use plasma losses. This legacy distinction is preserved.
- Table precision, catalogue-index-10 calorimetry scaling and every physical
  parameter remain as in Week 3. No calibration or scientific correction is added.

An independent distance for zero-redshift objects is not implemented; observer
conversion rejects z<=0 instead of generating infinities. This does not affect
the 11 current catalogue rows. CANDELS remains out of scope.

## Running and checking

In the astro environment:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
PYTHONPATH=src python -m congruents input/cat_nt.txt output/week4-smoke \
  --threads 4 --tables 16 --cosmic-rays 16 --photons 8 --cells 16
PYTHONPATH=src python -m congruents input/cat_nt.txt output/production \
  --threads 4 --cache data/python-cache
python tests/compare_outputs.py C_OUTPUT PYTHON_OUTPUT
```

The comparison command checks every standard output for existence, shape,
NaNs, infinities, exact reference zeros and numerical tolerance, returning
nonzero on failure. Default rtol=3e-6, atol=1e-280 allows for numerical integration
and six-digit output. This is not a percentage-similarity score and does not
exclude small values from the pass/fail comparison.

The independent reference drivers compile the original C source/headers,
changing only resolution and adding raw exports. The Week-4 reference continues
through the entire original output stage. It is not a wrapper around Python or
the new adapter. Tests use fresh temporary directories.

## Validation evidence and remaining gates

- The combined local suite passes all 41 tests (56.6 seconds in the recorded
  run). Both complete exported-file comparisons pass at the original tolerance.
- All 11 galaxies pass in-memory and exported-file comparisons at 16x16 tables /
  16 CR energies / 8 photons / 16 cells and at 32x32 / 32 / 16 / 32.
- Geometry and synthetic resampling match original C/GSL at rtol=1e-12 for all
  catalogue redshifts plus z=0.1, 1 and 3.
- All three bundled EBL tables match GSL at nodes and extrapolated queries,
  including below-grid energies, high redshift, energy capping and clipping.
- Tests cover explicit-tau validation, read-only/borrowed arrays, existing-output
  refusal, native failure recovery, radio ABI and one/four-thread equivalence.
- One 1000x1000 IC-3000 emission table generated successfully in Python in
  54.5 seconds locally. This alone does not validate the full production model.
- Full-size Python and independently executed original C catalogue runs both
  completed: 1000x1000 tables, 1000 CR energies, 500 photon energies, 500 cells,
  all 11 galaxies. All 39 standard files pass, covering 326,550 values.
  Only two printed values differ (both in CR_specs.txt); maximum relative
  error is 1.3698288624609998e-7, below the unchanged 3e-6 tolerance. See
  `WEEK4_PRODUCTION_COMPARISON.json` for the per-file report.
- The first Python run took 5560.48 seconds (92.67 minutes), generating 25
  tables and reusing the one previously timed IC-3000 table. It ran alongside
  the C reference for part of that period; this is not an isolated speed
  comparison.
- The cached repeat completed in 1955.95 seconds (32.60 minutes), with 26 cache
  hits and zero misses. All 39 files again pass the C comparison and are
  numerically identical to the first Python run at exported precision. These
  local wall times are observations, not a controlled speedup benchmark.
- Resident memory was observed at approximately 1 GB during the cached run;
  this is not a peak measurement. The macOS `/usr/bin/time -l` wrapper printed
  elapsed/user/system time but failed its `kern.clockrate` resource query with
  `Operation not permitted`, exiting 1 after the application completed and
  exported its results. Peak-memory statistics are therefore unavailable.
- Reproduction at the reference resolution is not a mathematical convergence
  proof or independent validation of the reference physics. No convergence
  threshold has been asserted from these agreement checks.
- A Linux full-suite CI job is configured alongside the existing portable
  macOS/Linux job, but has not been executed here. No Linux/HPC pass is claimed.
  GNU nested callback/trampoline restrictions must be checked on target HPC.
- The earlier focused runtime sanitizer/leak audit is not full-production
  sanitizer coverage. A clean scientific regression is not a memory proof.

Week 4 has not been pushed or merged. The Week-3 merge authorization does not
authorize a Week-4 merge.
