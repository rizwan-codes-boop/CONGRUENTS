# Week 1: C-library boundary and Python connection

Implementation date: 2026-09-06. Branch: feature/c-library-api.

## Scope

Week 1 is a working interface foundation, not an end-to-end Python galaxy
solver. The existing physics headers and production drivers are unchanged.
The pilot directly includes the existing ionisation header; it does not copy
the equation into Python or introduce new physical constants.

The supervisor's SLUG example is followed in loading a shared library,
declaring every ctypes restype/argtypes pair, and keeping C allocations behind
an opaque pointer. The C source has not been reorganised: moving the large
physics implementation now would unnecessarily broaden the regression scope.

## Daily deliverables

| Day | Delivered | Review evidence |
|---|---|---|
| 1 | Public ABI, units, buffer ownership, status codes | csrc/include/congruents.h |
| 2 | Shared-library target alongside original executables | make shared; make check |
| 3 | ctypes loader, ABI check, explicit signatures | src/congruents/_bindings.py |
| 4 | C-owned context, deterministic close, recoverable pilot errors | lifecycle, raw-boundary and error-recovery tests |
| 5 | Existing numerical function called through Python | comparison with independent direct-C executable, 1/4 threads |

## Run locally

From CONGRUENTS-pyc:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c check test-week1
PYTHONPATH=src python3 examples/week1.py
```

DEPENDENCY_ROOT points to the already installed compiler and libraries,
read-only. No dependency symlinks, tables or outputs are copied from the
baseline. Binaries and library products go into this checkout's bin/ and
build/, both ignored by Git.

Optional editable Python installation:

```sh
python3 -m pip install -e .
```

The Python packaging step does NOT compile or bundle the C library. Source/
editable use locates build/libcongruents.dylib (macOS) or .so (Linux). For a
separately installed Python package, set CONGRUENTS_LIBRARY to the absolute
library path or pass library= to Context. Distributable wheels are not a
Week-1 deliverable.

The tested toolchain is local GNU GCC 16.1 on macOS ARM64 with OpenMP.
The shared-link flag selects .so/-shared on Linux, but the inherited dependency
configuration is macOS-specific: Linux builds require toolchain/link-variable
overrides and are not yet validated. A fresh portable dependency installation
is still future work.

## Public API contract (ABI version 1)

| Function | Inputs | Output/ownership |
|---|---|---|
| cg_version | none | borrowed static UTF-8 version string |
| cg_abi_version | none | unsigned ABI version |
| cg_openmp_enabled | none | whether library was compiled with OpenMP |
| cg_status_message | integer status | borrowed static error message |
| cg_context_create | positive C-int thread count, output handle pointer | caller receives C-owned handle; status returned; failure sets handle NULL |
| cg_context_destroy | valid handle or NULL | releases context; non-NULL handle must be destroyed only once |
| cg_ionisation | borrowed context, count, double energy buffer, density, separate double output buffer | status; caller owns output; C retains neither buffer |

Status 0=success, 1=invalid argument, 2=allocation failure, 3=non-finite
numerical result. No library function changes GSL's global error handler.
The pilot does not call GSL and has no abort/exit error paths.

Energies are TOTAL electron energies in GeV, strictly above the existing C
electron rest-energy constant. Density is finite, non-negative hydrogen
number density in cm^-3. Loss is negative GeV/s; zero density produces zero.
Bounds validate the numerical API, not the physical accuracy of the empirical
loss formula throughout every conceivable energy range.

Zero-length calls are supported. Nonzero buffers must contain count doubles
and must not overlap. Raw C callers remain responsible for valid pointer
lifetimes; arbitrary dangling pointers cannot be made safe by ctypes.
Invalid-argument errors leave output untouched; discard the whole output for
any numerical-result error. Creation should receive a fresh output slot, not
overwrite an existing live handle.

Python copies iterable inputs into contiguous ctypes double buffers and returns
a list. This avoids a NumPy dependency for the pilot. It translates status codes
to exceptions, serialises calls/close on the same Python Context, and makes
close idempotent. Prefer "with Context(...)"; finalisation is only a fallback.
Distinct contexts may run concurrently. The native batch loop uses OpenMP
num_threads for that context, without mutating the global thread setting.

## Verification

- Both original production executables build and pass usage-message smoke tests.
- Direct-C versus ctypes comparison: 5 energies × 4 densities × 2 thread
  settings, relative tolerance 1e-14; zero results checked without a nonzero
  absolute tolerance.
- 1,000 create/use/free cycles; double Python close and use-after-close checks.
- Invalid energy, density, thread count, missing library and raw NULL argument tests.
- Numerical overflow returns an exception and a subsequent valid call succeeds.
- Concurrent independent contexts give identical results across thread counts.

These checks are not a memory-leak detector, a full-galaxy regression, a speedup
benchmark, or a proof of paper agreement. Existing production compiler warnings
(signedness, unused variables and possible uninitialised legacy arrays) remain.
The linker also warns about duplicate emutls linkage in the inherited toolchain.

## Next boundaries (not yet implemented)

Week 2 adds named galaxy/configuration inputs, units at the Python boundary,
table-context ownership and derived galaxy properties. Week 3 exposes the
solver/emission. Full-solver aborting error paths must be adapted before those
functions are callable from Python; Week-1 safety applies only to the exported
pilot, not arbitrary legacy functions. Astropy-based unit convenience can be
added without replacing the baseline C constants.

Physics corrections, including the issues recorded in the component map, must
remain separate from interface changes. Do not regenerate baseline tables or
silently change equations while extending the wrapper.

No commits, pushes, tags, pull requests or merges are performed by this task.
