"""Reduced-grid integration check, not a scientifically converged spectrum."""
from pathlib import Path
from congruents import Preparation, Grid, SolverGrid, load_catalogue, solve

root = Path(__file__).resolve().parents[1]
with Preparation(load_catalogue(root/"input/cat_nt.txt"), Grid(16,16)) as preparation:
    result = solve(preparation, SolverGrid(cosmic_rays=16, photons=8, cells=16), threads=4)
    print("Galaxies:", len(preparation.galaxies))
    print("Electron components:", list(result.electrons))
    print("Emission components:", list(result.emission))
    print("IC primary disc shape:", result.emission["IC_primary_disc"].shape)
    print("Reduced grid only: not production results.")
    # Optional explicit export: result.save("week3-results.npz")
