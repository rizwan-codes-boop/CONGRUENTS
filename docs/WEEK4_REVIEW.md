# Week 4: observer workflow and production acceptance

Status: started, not complete. Branch: `feature/week4-observer-workflow`.
Week 3 was published and merged into main at `73d64d5` on user request.

## First increment: implemented

`src/congruents/observer.py` contains Python-only luminosity distances, the
reference distance factor, and one-component observer SED conversion. There
are no new C production kernels or callbacks. The standalone C test executable
is only an independent reference check.

- Cosmology follows the actual numbers in `cosmo_params.h`: H0=70 km/s/Mpc,
  Om=0.3, Ode=0.7, flat, without radiation. Its Planck-2018 comment does not
  match those numbers; this port preserves the numbers, not that label.
- Distance follows `cosmo_funcs.h:d_l_MPc`, with Python quadrature. The
  `distmod` factor in `spectra.c:1003` is cm^-2, not a magnitude modulus.
- Resampling follows `spectra.c:1028`: linear interpolation on E/(1+z), onto
  the original photon grid. Values outside the contracted domain are zero.
- Output follows `file_io.h:write_2D_spec_file`: multiply by exp(-tau_internal),
  exp(-tau_EBL), distance factor and E². Output is GeV cm^-2 s^-1.
- Both optical-depth arrays are required. No default zero attenuation and no
  silent broadcasting. Explicit zeros are allowed only as an intentional
  unattenuated diagnostic. Inputs remain unchanged; results are read-only.

The optical-depth convention intentionally follows the existing writer's
array indexing. This is not a claim that indexing is physically preferable.
Any correction to redshift/attenuation evaluation belongs in a separate change.
The generic transform does not select which particle components should receive
attenuation; production component selection is not connected yet.

## Validation

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test-week4 PYTHON=python
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
```

The initial tests compare distances and a synthetic component transform against
the original C cosmology and GSL linear interpolation for all 11 catalogue
redshifts and z=0.1, 1 and 3, at rtol=1e-12. They also check the zero endpoint,
required attenuation, invalid inputs and explicit unattenuated diagnostics.
This is not a comparison of final galaxy spectra or production-resolution runs.
Zero-redshift flux conversion is rejected: an independent-distance API remains
future work rather than returning infinities.

Local macOS validation: all three initial Week-4 tests and all 36 tests in the
combined regression suite passed. This does not establish Linux/HPC coverage.

## Remaining increments

1. Port and independently verify the actual EBL-table reader/interpolation,
   including extrapolation, energy cap, nonnegative clipping and provenance.
   Implement internal gamma-gamma optical depth and verify all galaxies.
2. Connect source results to all observer components; add totals, legacy-format
   export, spectral integration and the missing loss/radio/energy diagnostics.
   Compare every applicable output against a full independent C driver.
3. Run full-resolution tables and 500-cell production comparisons, convergence,
   memory and performance checks; exercise Linux/HPC builds and document limits.

Keep the Week-3 Python-prepared transport/free-free boundary. Do not put serial
postprocessing or table generation back in C. No CANDELS integration or physical
parameter corrections are included. This Week-4 branch is not authorized for
merge by the user's Week-3 merge request.
