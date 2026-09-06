# Reconstruction record

## Scope

This directory is a C-only reconstruction of the released CONGRuENTS source.
The production physics headers and `spectra.c` remain in C. OpenMP continues
to parallelise the galaxy and interpolation loops in compiled code.

The immediate target is reproducibility of the MNRAS 2023 calculation. It is
not yet a revised-physics fork: known scientific discrepancies should first be
measured against the archived paper products before any formula is changed.

## Files missing from the source snapshot

The snapshot referenced, but did not contain, these project headers:

- `physical_constants.h`
- `astro_const.h`
- `cosmo_params.h`
- `math_funcs.h`
- `gsl_decs.h`

They have been reconstructed from call sites and unit suffixes. Compatibility
forwarders for the source's non-standard flat GSL includes are also present.
`omp.h` supplies a serial fallback, while the normal Makefile build selects the
real OpenMP runtime.

## Explicit reconstruction assumptions

- Physical constants use CODATA 2018 and IAU nominal conversion values.
- Cosmology is flat Planck 2018: H0 = 67.4 km/s/Mpc, Omega_m = 0.315,
  Omega_Lambda = 0.685, Omega_k = 0.
- The reconstructed 2-D spline wrapper uses GSL bilinear interpolation. This
  matches the explicit `gsl_interp2d_bilinear` choice in `spectra.c`.
- The original GSL requirement was 2.7.1; the local build is GSL 2.8. A
  numerical comparison is still required to quantify any version sensitivity.

## Dependencies recovered locally

- GSL 2.8 is built statically under `.deps/`.
- Steven G. Johnson's cubature source is under `vendor/cubature`, with a static
  library under `.deps/`.
- GNU GCC 16 and the OpenMP runtime are stored project-locally because the
  released code uses GNU nested functions unsupported by Apple Clang.

## Known legacy issues deliberately not corrected yet

These are recorded so a successful legacy reproduction is not confused with a
paper-faithful corrected model:

- `CR_spectra/gal_rad.h`: the dilution coefficient contains `4*pi`, whereas
  the corresponding paper expression uses `2*pi`.
- `data_calc.h`: dust-temperature table interpolation applies the fractional
  weights in the opposite order from conventional lower/upper interpolation.
- `data_calc.h`: a halo diagnostic labelled for plasma losses evaluates the
  ionisation lifetime.
- `spectra.c`: index-specific corrections involving galaxy index 10 are
  hard-coded and require provenance from the original production catalogue.
- `data_objects.h`: floating-point comparisons use integer `abs` rather than
  `fabs`; this can make cache validation accept incompatible interpolation
  tables. It is preserved for the first legacy-output comparison.

## Reproduction stages

1. Build and smoke-test both executables with `make check`.
2. Generate expensive lookup tables with `make run-precompute`.
3. Produce the complete spectra with `make run`.
4. Compare every output table against the archived paper-production outputs.
5. Only after the legacy baseline matches, apply scientifically corrected
   formulas one at a time and report their effect.

Interpolation generation is intentionally a separate stage because the
1000-by-1000 radiative tables can take substantial time even with OpenMP.

## Runtime corrections applied during validation

- GSL interpolation accelerators associated with the shared bremsstrahlung
  and synchrotron tables are now private to each OpenMP galaxy iteration. The
  accelerator is a mutable cache and sharing it caused data races.
- `q_pi` now returns zero when the proton energy lies outside the calorimetry
  table. The legacy `q_nu` integration bounds otherwise request values as low
  as approximately 1e-15 GeV on a table beginning at 1e-3 GeV, causing GSL to
  abort. This guard enforces the table's physical domain without extrapolation.
- The redshifted final-spectrum interpolation now returns zero above its
  contracted tabulated energy domain. The legacy loop evaluated every point
  up to 1e8 GeV even though each redshifted spline ends below 1e8 GeV.
- The production target creates `output/tau_loss/`, which the source assumes
  exists when writing the loss-time diagnostics.
