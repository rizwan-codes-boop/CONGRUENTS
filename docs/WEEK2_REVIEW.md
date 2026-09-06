# Week 2: inputs, galaxy properties and tables

Implementation date: 2026-09-06.
Branch: feature/model-inputs-and-tables.
Week 1 was merged and pushed to main as 4ce7c98 before starting this work.

## What is delivered

| Day | Feature | Main files / tests |
|---|---|---|
| 1 | Named, validated four-column catalogue; round-trip writer; optional Astropy quantity conversion | inputs.py; catalogue and Astropy tests |
| 2 | Native galaxy-property batch function, OpenMP in C | csrc/preparation.c; all 11 gal_data rows compared with saved baseline |
| 3 | Python-owned NumPy IC/BS/SY arrays, native generation/evaluation, Python cache | preparation.py and C adapter; all 14 families and cache tests |
| 4 | Per-galaxy photon fields, integrated densities, combined emission/Gamma IC tables | direct legacy combination and Urad_Ub comparisons |
| 5 | End-to-end preparation example, direct-C reference executable, regression fixtures | examples/week2.py; tests/direct_preparation.c; tests/test_week2.py |

No electron solve or final galaxy spectrum is exposed yet. No physical equations,
constants, production drivers, original table routines or FF/IC corrections have
been modified. The interface calls existing helpers. Galaxy-property expressions
that were embedded inside the driver are mirrored in the adapter, retaining
the original constants and operation order, with regression tests guarding drift.

## Run

From CONGRUENTS-pyc, in the existing astro environment:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c check test-week1 test-week2
PYTHONPATH=src python examples/week2.py
PYTHONPATH=src python examples/week2.py
```

If make selects a different Python interpreter:

```sh
make DEPENDENCY_ROOT=../CONGRUENTS-c PYTHON=/Users/rizzy/anaconda3/envs/astro/bin/python check test-week1 test-week2
PYTHONPATH=src /Users/rizzy/anaconda3/envs/astro/bin/python examples/week2.py
```

The example uses 8×8 tables: quick regression checks, NOT production-accuracy
physics. On a cold cache it generates 26 planes; on the second run all 26 load
from cache. The current catalogue needs 8 fixed-field IC planes, 4 CMB planes,
12 FIR planes, one BS plane and one SY plane. This is 14 logical families.

Grid() defaults to the unchanged 1000×1000 plane resolution. Increasing grids
can be expensive: do not use the smoke-example values to interpret the paper.
Only four selected full-resolution legacy IC cells were reintegrated during
this task, not the full production grids or spectra.

## Python API

```python
from congruents import Galaxy, Grid, Preparation, load_catalogue

galaxies = load_catalogue("input/cat_nt.txt")
with Preparation(galaxies, grid=Grid(8, 8), threads=4,
                 cache="data/my-preparation") as prep:
    properties = prep.properties
    radiation = prep.radiation(0, [1e-12, 1e-9])
    prep.prepare_tables()
    with prep.combined_ic(0, "emission") as kernel:
        rates = kernel.evaluate([1e-6], [1.])
```

Python owns the persistent NumPy table buffers; C borrows their addresses for
numerical evaluation. Preparation owns its field tables and closes their C
descriptors. table() returns one of these
owned tables; callers should generally not close it independently before reusing
the preparation. combined_ic() returns a NEW independent handle which callers
must close. Snapshots are independent Python copies. Closing is idempotent;
use-after-close raises an exception. Methods serialize calls/close per owner.
The caller must not mutate preparation properties/configuration while in use.

Table.x, Table.y and Table.values expose read-only NumPy views; values has
shape (ny, nx), dtype float64 and C-contiguous layout. Use values.copy() for
editable analysis data. Strong references keep buffers alive for the entire C
descriptor lifetime. C never frees Python-owned buffers; retained NumPy views
remain valid after the descriptor is closed. Do not bypass the read-only flags
or mutate private buffers while C is using them.

Cold generation and combination still use the unchanged C builders' temporary
arrays. Table construction copies those results once into NumPy and immediately
releases the temporary native table. Warm-cache loading goes straight to NumPy.
The only persistent table-side C allocation is a small borrowed-array descriptor.
This avoids rewriting legacy kernels merely to eliminate their scratch buffers.

The optional Astropy extra only handles units, not replacement physics constants:

```python
from astropy import units as u
g = Galaxy.from_quantities(
    redshift=0.001, stellar_mass=1e10*u.Msun,
    radius=2000*u.pc, sfr=1*u.Msun/u.yr)
```

Install the optional dependency with pip install -e '.[units]' if required.
The supplied astro environment already provided Astropy 7.2.0 and passed this test.
NumPy is now required. Astropy remains optional; its unit test skips if absent.

## Units and layout

- Input catalogue has FOUR columns: redshift, stellar mass [Msun], effective
  radius [kpc], SFR [Msun/year]. Preserve row order and duplicate model cases.
  This corrects the earlier component-map wording that said five columns.
- Properties are ten named values: the nine gal_data columns plus halo B.
- Photon fields: CMB, FIR, 3000, 4000, 7500, UV, total, in cm^-3 GeV^-1.
- Integrated densities: total, magnetic, CMB, FIR, 3000, 4000, 7500, UV,
  in eV/cm³, matching Urad_Ub ordering.
- IC emission: x=emitted photon GeV, y=TOTAL electron GeV, z=1/(s GeV).
- IC Gamma: x=electron energy transfer DeltaE GeV, y=TOTAL electron GeV,
  z=1/(s GeV). This is not itself an energy-loss rate.
- BS: x=photon GeV, y=TOTAL electron GeV, z=mb/GeV.
- SY: x dimensionless, y is a one-element unused sentinel, z dimensionless.
- Table snapshots flatten z as iy*nx+ix; no implicit transpose.
- Evaluation is linear/bilinear in native coordinates, zero outside the table.
  It does not interpolate logarithms of rates.

## C API and safety

The additions retain ABI 1 and do not alter Week-1 signatures:
cg_galaxies, cg_preparation_bounds, cg_radiation, cg_table_create,
cg_table_import, cg_table_borrow, cg_table_destroy, cg_table_shape, cg_table_copy,
cg_table_eval, cg_temperature_grid and cg_combine_ic. All signatures are
explicit in the public header and preparation.bind().

cg_table_borrow creates a descriptor without copying data. Its caller must
provide valid contiguous double buffers of lengths nx, ny, nx*ny and retain
them until destruction. cg_table_import remains available for C-owned scratch
copies, but normal Python cache loading does not use it.

C validates finite inputs, energy bounds, shape limits, monotonic axes and
finite/non-negative imported rates. Negative/zero masses, radii or SFRs are
rejected because the current empirical pipeline uses logarithms. The exposed
preparation domain restricts redshift to 0..20 and generated tables to the
documented numerical ranges in cg_table_create; these are interface guards,
not newly inferred physical validity ranges.

The adapter's own allocations are checked and return status codes. However,
legacy builders and GSL wrappers still contain unchecked internal allocations
and process-aborting paths. Valid small-grid runs and rejected boundary inputs
are tested; this is NOT a guarantee against a process failure under memory
exhaustion or every extreme legacy numerical failure. No global GSL error
handler is changed. Broader hardening remains required before claiming a fully
fault-contained production shared library.

## Cache contract

New .cgt files are deliberately separate from legacy .txt caches. They store:
magic/schema, JSON metadata, dimensions, little-endian double arrays and SHA-256
payload checksum. Metadata includes the loaded library's checksum, grids,
temperature, field, table type, axes and units. Changes to the binary/configuration
produce new keys. Writes are atomic; damaged or mismatched caches fail explicitly.
No pickle or executable deserialisation is used.

Python owns table arrays and handles cache I/O; C computes the numerical values.
Original text caches are NOT imported by this API. Regeneration or a future
explicitly validated importer is required to use production tables here.
Caches are local build artifacts, not intended as cross-toolchain reproducibility
proof. The library hash is conservative: even an unrelated rebuild can invalidate
reuse. External linked runtime versions are not fully fingerprinted yet.

## Preserved behaviours and explicit limitations

1. The existing reversed temperature weights are retained:
   lower*frac + upper*(1-frac). A test explicitly checks that a lower-node
   fraction of zero selects the upper plane. This is baseline fidelity, not
   a claim that the interpolation is physically correct.
2. Temperature grids retain the legacy count/rounded-bound policy. Fewer than
   two planes (e.g. a single identical dust temperature) are rejected instead
   of entering the legacy unsigned-loop failure. The upper temperature
   endpoint, which has no next plane, is also rejected. Use a catalogue with
   adequate coverage; automatic padding would be a separate policy change.
3. Individual 3000/4000 radiation diagnostic helpers use their legacy luminosity
   choices. Do not assume the separately reported components sum to total;
   regression follows the existing output rather than repairing it.
4. Galaxy-property batches use OpenMP in C. Table construction remains serial
   as in the inspected legacy loops; no Python multiprocessing is added.
5. No CANDELS integration, diffusion-grid/injection solve, EBL application or
   final photon spectra are introduced in Week 2.
6. Compiler warnings in unchanged source remain; no full numerical convergence
   study, memory-leak audit or fresh-machine packaging validation is claimed.

## Verification summary

On macOS ARM64 / GCC 16.1 / astro Python:

- Both existing production executable smoke checks pass.
- All seven Week-1 tests pass.
- All twelve Week-2 tests pass, including Astropy (not skipped).
- NumPy ownership/read-only views, retained-buffer lifetime, and zero-copy
  native borrowing pass explicit tests.
- All 11 galaxy-property and radiation-density rows agree within 6e-7,
  accounting for saved-file precision.
- All 26 small-grid planes match the independently compiled legacy builders.
- Combined emission/Gamma IC matches the original construct_IC_gso2D for four
  catalogue cases within 2e-13.
- Cold/warm cache payloads are exactly equal; corruption is rejected.
- Four original 1000×1000 IC cells agree with fresh integration within 1e-5,
  allowing for rounded stored energy coordinates.

Review the code and these tests on feature/model-inputs-and-tables before merging
Week 2. This feature branch is kept separate from main.
