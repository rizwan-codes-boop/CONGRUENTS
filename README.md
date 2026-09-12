# CONGRUENTS — Python/C development and preserved C reference

The interface uses Python for serial preparation and C/OpenMP for galaxy batches,
following the supervisor-confirmed [requirements](docs/PROJECT_REQUIREMENTS.md).
The existing Week-1/Week-2 features have been revised. Week 3 exposes the
two-zone solver and source emission; observer-frame production is still pending.

See [Week 3: solver, emission and tests](docs/WEEK3_REVIEW.md). Its optional
native library requires the existing GSL/cubature dependencies:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c test-week3 PYTHON=python
PYTHONPATH=src python examples/week3.py
```

This is a reduced-grid demonstration for all 11 galaxies, not a converged
production spectrum. No final Figure 9 or observer-frame agreement is claimed.

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

The solver uses ABI 2: rebuild with `make solver` after updating. See
[`docs/WEEK3_REVIEW.md`](docs/WEEK3_REVIEW.md) for the remaining native helpers
and the limits of this precomputation boundary.

## Building on other systems

Week 4 has started with Python observer-frame distance and resampling primitives.
Attenuation generation and final production outputs are not connected yet; see
[`docs/WEEK4_REVIEW.md`](docs/WEEK4_REVIEW.md) for scope and tests.

The interface library requires GNU GCC with OpenMP, but **not** GSL/cubature:

```sh
make shared test-portable CC=gcc PYTHON=python
```

On macOS with Homebrew GCC, pass the versioned compiler path, for example
`CC="$(brew --prefix gcc@14)/bin/gcc-14"` after installing that formula.
On HPC, load the site's GCC/Python modules and pass its compiler path. Restrict
native threads to allocated CPUs. macOS/Linux CI is configured for the portable
subset; successful Linux/HPC execution has not yet been established locally.

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
