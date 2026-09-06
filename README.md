# CONGRuENTS — reconstructed C production tree

## Inputs and table preparation (Week 2)

See [docs/WEEK2_REVIEW.md](docs/WEEK2_REVIEW.md) for the catalogue, optional
Astropy units, C galaxy properties/radiation and Python-owned NumPy IC/BS/SY tables.
Python handles storage, caching and orchestration; C performs numerical work.
Run the small-grid example twice to verify cache reuse:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c check test-week1 test-week2
PYTHONPATH=src python examples/week2.py
```

Use the astro Python environment for the optional unit tests. The example is
preparation-only at smoke-test resolution, not a final galaxy spectrum.

## Python–C development (Week 1)

The pilot shared library and Python bindings are documented in
[docs/WEEK1_REVIEW.md](docs/WEEK1_REVIEW.md). From this development checkout,
with the existing sibling dependency installation:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c check test-week1
PYTHONPATH=src python3 examples/week1.py
```

This exposes the unchanged C ionisation function, not yet the complete galaxy
solver. Physics remains in C, including the pilot's OpenMP batch loop.

CONGRuENTS produces cosmic-ray, neutrino, gamma-ray, and radio non-thermal
spectra for galaxy catalogues. This folder preserves the end-to-end model in C
and restores the files and build dependencies absent from the available source
snapshot.

## Quick start

From this directory:

```sh
make check
make run-precompute
make run
```

`make check` compiles the two production executables and performs lightweight
command-line smoke tests. `make run-precompute` generates the expensive inverse
Compton, bremsstrahlung, and synchrotron lookup tables in `data/`. `make run`
reads those tables and writes the complete result set to `output/`.

Set the C parallelism with, for example:

```sh
OMP_NUM_THREADS=8 make run-precompute
OMP_NUM_THREADS=8 make run
```

## Executables

- `bin/create_interp_objects input/cat_nt.txt data` builds reusable numerical
  interpolation tables.
- `bin/spectra input/cat_nt.txt data output` runs the full production model.

Both are compiled with GNU C and OpenMP. No model component has been converted
to Python.

## Directory layout

- `CR_spectra/` — inverse Compton, synchrotron, bremsstrahlung, ionisation,
  diffusion, and radiation-field physics.
- `input/` — recovered galaxy catalogue and EBL optical-depth tables.
- `data/` — generated interpolation tables; initially empty.
- `output/` — generated production spectra; initially empty.
- `vendor/` and `.deps/` — recovered source dependencies and local libraries.
- `RECONSTRUCTION.md` — missing-file provenance, assumptions, known legacy
  issues, and the staged validation plan.

## Dependency reconstruction

The original repository named GSL and cubature but supplied neither a build
system nor its five private headers. This copy includes locally built GSL 2.8,
vendored cubature, and a project-local GNU compiler/OpenMP runtime. GNU C is
required because the original implementation defines numerical callbacks as
nested functions, a GNU extension rejected by Apple Clang.

The project-local compiler packages are intentionally kept out of the physics
source. The Makefile discovers them beneath `.conda-pkgs/`.

## Scientific status

Compilation is only the first reproducibility gate. The restored constants,
cosmology, and spline wrappers are documented assumptions until numerical
outputs are compared against the original paper-production files. See
`RECONSTRUCTION.md` before interpreting results or changing formulas.
