# Paper-faithful implementation fork

`CONGRUENTS-i` is an experimental fork of `CONGRUENTS-c`. Its purpose is to
measure the output changes produced by making the implementation agree with
the equations and stated corrections in Roth et al. (2023), MNRAS 523,
2608-2629.

## Applied corrections

1. Radiation-field dilution uses `2*pi` rather than the legacy `4*pi`
   denominator (Section 2.3).
2. The free-free Gaunt-factor branch changes at 1 GHz (equation 45).
3. The extra division by `f_e` is removed after the `f_e n_i^2` cancellation
   in the free-free optical depth (equations 44 and 51).
4. The special corrected-LMC model retains the paper's proton-calorimetry
   adjustment but no longer reduces primary-electron injection by a factor of
   ten (Section 5.2.2).
5. Lookup-table validation uses floating-point `fabs`, and the 3-D field upper
   bound is checked with a tight tolerance rather than the legacy `1e6`.
6. Temperature-dependent CMB/FIR lookup generation selects the kernel using
   the requested table type (`j`) rather than the radiation-field index (`i`).
   The legacy index mix-up could write an emission table with the transition-
   rate kernel, or vice versa.

The expensive IC lookup tables must be generated inside this directory. The
local compiler and dependency trees are shared read-only with `CONGRUENTS-c`;
model source, lookup tables, binaries, and outputs remain separate.
