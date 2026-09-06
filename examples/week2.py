"""Small-grid integration example; these tables are NOT production resolution."""
from pathlib import Path
from congruents import Grid, Preparation, load_catalogue

root = Path(__file__).resolve().parents[1]
galaxies = load_catalogue(root / "input/cat_nt.txt")
with Preparation(galaxies, grid=Grid(8, 8), threads=4,
                 cache=root / "data/pyc-week2-example") as model:
    print("Galaxies:", len(galaxies))
    print("First galaxy:", model.properties[0])
    print("Radiation densities [eV/cm³]:",
          model.radiation(0, [1e-12, 1e-9])["urad_ub_ev_cm3"])
    print("Prepared planes:", model.prepare_tables())
    print("Cache hits:", model.cache_hits, "generated:", model.cache_misses)
    with model.combined_ic(0, "emission") as emission:
        print("Python-owned NumPy table:", emission.values.shape,
              emission.values.dtype, "writable:", emission.values.flags.writeable)
        print("IC rate [1/(s GeV)] at E_gamma=1e-6 GeV, E_e=1 GeV:",
              emission.evaluate([1e-6], [1.])[0])
    with model.combined_ic(0, "gamma") as transition:
        print("IC transition rate [1/(s GeV)] at DeltaE=1e-6 GeV, E_e=1 GeV:",
              transition.evaluate([1e-6], [1.])[0])
    print("Smoke-test grid: 8×8; do not use these values as converged spectra.")
