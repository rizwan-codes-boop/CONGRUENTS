# Standing project requirements

Confirmed by the user following supervisor discussion, 2026-09-11.
These supersede the earlier proposal to retain a broad C numerical backend.
Apply throughout the project unless the user explicitly changes them.

## Supervisor decisions

1. **Python–C split:** Everything outside the OpenMP loops belongs in Python.
   Prepare the native C variables/buffers before entering the OpenMP loops to
   avoid repeated conversion and Python–C crossing overhead. Do not interpret
   "computationally expensive" alone as permission to retain a component in C.
2. **Reference implementation:** Use CONGRUENTS-i, not the original incomplete
   astromatt42 checkout or a different experimental fork.
3. **Parameters and behaviour:** Keep these the same as CONGRUENTS-i. Do not
   introduce new scientific parameter choices or model behaviour as part of
   the interface conversion.
4. **Initial task:** Reproduce the same task/workflow as CONGRUENTS-i. New use
   cases, CANDELS integration and model extensions are deferred until requested.
5. **Validation:** Check all reference cases and applicable outputs/components,
   not only selected galaxies, ionisation or a few lookup-table cells. Existing
   smoke tests and spot checks are useful but not final acceptance.
6. **Parallelism:** Parallelise over galaxies only. Do not add parallel
   lookup-table generation, energy-bin processing or parameter scans.
7. **Platforms:** Aim for macOS, Linux and HPC environments. Local Mac success
   is not evidence of portability; document and test platform coverage.

Items 3–7 interpret the numbered answers against the preceding questions about
configurability, use case, validation, parallelism and target platforms.

## Implementation consequences

- Python owns inputs, units, configuration, data storage, caching, orchestration
  and output/plot handling, plus calculations outside the OpenMP regions.
- Native inputs must be prepared before parallel execution, with explicit
  dtype, shape, layout, units and lifetime contracts. Workers must not repeatedly
  call Python for conversions or table lookups.
- OpenMP loop bodies necessarily need executable native work. The exact
  placement of helper functions called by those bodies must be made explicit
  in the next architecture review; do not silently use that fact to retain
  the entire existing C backend.
- Inventory each routine against the active CONGRUENTS-i OpenMP call graph.
  In particular, serial lookup-table generation cannot remain in C merely
  because it is expensive under the revised requirement.
- Preserve CONGRUENTS-i physics, constants, normalisations, units, catalogue
  order and documented implementation behaviours. Propose scientific fixes
  separately; do not silently apply them during language conversion.
- Compare all 11 current catalogue cases and applicable production outputs.
  Agree and document numerical tolerances; exact bitwise equality is not
  assumed when numerical implementations change.
- Retain the small tests, but add end-to-end production validation before
  claiming completion of the interface.

## Current work and transition

At the time of the decision, Week 1 was merged and Week 2 was on
feature/model-inputs-and-tables, with serial numerical routines still in C.
The user subsequently authorised updating that work on a new branch and merging
it into main. The revision on refactor/python-serial-galaxy-openmp moves those
serial components to Python; see PYTHON_SERIAL_REVIEW.md for the reviewed
boundary, regression coverage and remaining full-solver acceptance gates.

Review the architecture against these decisions before extending the solver
interface. This document records requirements only: it does not authorise a
merge, push, scientific correction, or immediate source rewrite.

Keep existing user edits, especially the uncommitted loader-docstring change,
separate from implementation commits unless the user requests otherwise.
