# CONGRUENTS Python–C project instructions

Read docs/PROJECT_REQUIREMENTS.md before planning, editing, or reviewing this
project. Those supervisor-confirmed requirements supersede earlier architecture
plans and Week-1/Week-2 documentation where they conflict.

In particular: use CONGRUENTS-i as the unchanged scientific reference; put
everything outside the OpenMP loops in Python; prepare native variables before
parallel execution; parallelise galaxies only; validate all reference cases
and applicable outputs; target macOS, Linux and HPC.

The exact native loop-body/helper boundary needs an explicit architecture review.
Do not substitute "keep all expensive calculations in C" for the approved split.
The revised preparation boundary is documented in docs/PYTHON_SERIAL_REVIEW.md.
The optional solver/emission boundary is documented in docs/WEEK3_REVIEW.md.
Its ABI-2 revision precomputes transport and free-free arrays in Python; do not
move those calculations back into workers. Remaining native integration/solver
helpers are explicitly inventoried there, not assumed to be irreducible.
Full observer-frame orchestration and production acceptance remain outstanding.

Do not merge or push without user authorisation. Preserve unrelated user edits.
