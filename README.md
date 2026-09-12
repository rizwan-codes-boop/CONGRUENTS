# CONGRUENTS — Python/C development and preserved C reference

The interface uses Python for serial preparation and C/OpenMP for galaxy batches,
following the supervisor-confirmed [requirements](docs/PROJECT_REQUIREMENTS.md).
The existing Week-1/Week-2 features have been revised. Week 3 exposes the
two-zone solver and source emission. Week 4 adds Python attenuation, observer
components, diagnostics and legacy-format output; see the validation status below.

See [Week 3: solver, emission and tests](docs/WEEK3_REVIEW.md). Its optional
native library requires the existing GSL/cubature dependencies:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test-week3 PYTHON=python
PYTHONPATH=src python examples/week3.py
```

This example is a reduced-grid demonstration for all 11 galaxies, not a
converged production spectrum. Full-resolution observer outputs have separately
passed comparison against CONGRUENTS-i; see the Week-4 review. That is not an
independent comparison against digitized curves from the published paper.

See [the current architecture and validation review](docs/PYTHON_SERIAL_REVIEW.md).
The older Week-1/Week-2 reviews describe historical implementations.

## Interface quick start

In the existing astro environment, using the sibling compiler dependencies:

```sh
python -m pip install -e '.[units]'
make DEPENDENCY_ROOT=../CONGRUENTS-c test PYTHON=python
PYTHONPATH=src python examples/week1.py
PYTHONPATH=src python examples/week2.py
```

This builds library version 0.2.0, ABI 2. Rebuild it when upgrading from Week 2.
The small-grid example generates preparation tables and diagnostics, not final
galaxy spectra. Run it twice to see cold/warm cache behaviour.

Python owns ionisation diagnostics, grids, radiation fields, IC/BS/SY table
generation, interpolation, combination, storage and caching. The native library
contains the galaxy-property OpenMP loop, its directly used helpers, and minimal
context/ABI plumbing. Python also prepares normalization, calorimetry, diffusion,
primary/proton spectra and free-free arrays. The optional solver library adds
secondary-injection integrals, two-zone solves and nonthermal emission inside
galaxy loops. There are no Python callbacks
inside native workers.

The solver uses ABI 3: rebuild with `make solver` after updating. See
[`docs/WEEK3_REVIEW.md`](docs/WEEK3_REVIEW.md) for the remaining native helpers
and the limits of this precomputation boundary.

## Catalogue to observer outputs

```sh
# Small regression run, not a converged scientific prediction:
PYTHONPATH=src python -m congruents input/cat_nt.txt output/week4-smoke \
  --threads 4 --tables 16 --cosmic-rays 16 --photons 8 --cells 16
# Production counts (1000x1000 tables, 1000 CR bins, 500 photon bins, 500 cells):
PYTHONPATH=src python -m congruents input/cat_nt.txt output/production \
  --threads 4 --cache data/python-cache
```

Use a new output directory for each run. Existing outputs are never overwritten.
The production run can take substantially longer than the small tests. See
[`docs/WEEK4_REVIEW.md`](docs/WEEK4_REVIEW.md) for validation evidence and limits.

## Building on other systems

The interface library requires GNU GCC with OpenMP, but **not** GSL/cubature:

```sh
make shared test-portable CC=gcc PYTHON=python
```

On macOS with Homebrew GCC, pass the versioned compiler path, for example
`CC="$(brew --prefix gcc@14)/bin/gcc-14"` after installing that formula.
On HPC, load the site's GCC/Python modules and pass its compiler path. Restrict
native threads to allocated CPUs. macOS/Linux CI covers the portable subset;
a separate Linux full-suite job is configured but not yet verified here.
Successful Linux/HPC execution must not be inferred from local Mac checks.

The full regression suite also needs the preserved C reference executables and
their GSL/cubature dependencies. The portable subset is not a substitute for it.

## Preserved C production workflow

The files below remain the CONGRUENTS-i scientific reference. Their existing
standalone workflow is separate from the unfinished Python solver interface:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c check
make DEPENDENCY_ROOT=../CONGRUENTS-c run-precompute
make DEPENDENCY_ROOT=../CONGRUENTS-c run
```

The first command compiles and smoke-checks the C executables.
Precomputation writes legacy tables in `data/`; spectra writes into `output/`.
OpenMP parallelises galaxy processing. Lookup-table generation is serial.

- `CR_spectra/`, root C sources/headers: preserved reference physics/drivers.
- `src/congruents/`: Python interface and serial calculations.
- `csrc/`: narrow native galaxy-batch interface.
- `input/`: reference catalogue and EBL data.
- `tests/`: independent C comparisons, fixtures and portable checks.
- `RECONSTRUCTION.md`: missing-file provenance and reconstruction assumptions.

Generated caches, shared libraries and outputs are not committed. The interface
does not integrate CANDELS or change the reference scientific parameters.
