# Week-2 regression fixtures

These are small copies/excerpts of local CONGRUENTS-i generated outputs,
captured on 2026-09-06. They are NOT observations, digitised paper curves or
independent proof of the model's physics.

The two text files contain the header plus all 11 catalogue rows:

| Original relative to Code/ | SHA-256 |
|---|---|
| CONGRUENTS-i/output/gal_data.txt | 11c524cb61bf99c46b8a3948503d8190e509d40c6d6c1799b6651b9c8826c5f6 |
| CONGRUENTS-i/output/Urad_Ub.txt | 4576697395c1f8dd1521bc0808f5a9251dfebd0d9bfea6143fae243864ad1549 |
| CONGRUENTS-pyc/input/cat_nt.txt | 3a992ec673ae28d0089c2f1c318f8544c6f76fb5e00373581b4038b5b3f2ad8d |

production_ic_spot.json records four cells from a 1000×1000 IC_3000 emission
table: zero-based x indices 400/401, y indices 200/201, flattened y-major.
Its source checksum and saved (rounded) coordinates are embedded in the JSON.
The spot test reintegrates at those coordinates, not a coarse-grid interpolation.

Text outputs have seven significant digits. Property/radiation comparisons
use 6e-7 relative tolerance. The spot comparison allows 1e-5 for rounding of
both coordinates and values. Full small-grid direct-C comparisons use exact
table equality, with 2e-13 relative tolerance for combined IC round-off.

Historical output provenance before this snapshot is not independently
authenticated. For this reason the tests ALSO compile a separate executable
using the unchanged legacy table builders and combination function.
