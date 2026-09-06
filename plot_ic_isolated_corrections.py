"""Compare the two inverse-Compton corrections independently.

Ratios are measured against the unmodified CONGRUENTS-c reconstruction.
The combined curve comes from CONGRUENTS-i, which contains both IC fixes.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CODE = Path(__file__).resolve().parent.parent
ROOTS = {
    "Legacy": CODE / "CONGRUENTS-c" / "output",
    r"$2\pi$ only": CODE / "CONGRUENTS-2pi-only" / "output",
    "Lookup only": CODE / "CONGRUENTS-lookup-only" / "output",
    "Full corrected fork": CODE / "CONGRUENTS-i" / "output",
}
DEST = CODE / "CONGRUENTS-i" / "output"

GALAXIES = [
    ("NGC 2146", 4),
    ("NGC 2403", 3),
    ("M82", 5),
    ("NGC 253", 2),
    ("SMC", 9),
    ("M33", 7),
]
IC_FILES = (
    "spec_IC_1_z1.txt",
    "spec_IC_2_z1.txt",
    "spec_IC_1_z2.txt",
    "spec_IC_2_z2.txt",
)
ONE_KEV_GEV = 1.0e-6
X_LIMITS = (1.0e-16, 1.0e5)


def load(root, filename):
    return np.loadtxt(root / filename, skiprows=1)


def load_ic(root):
    return sum(load(root, filename) for filename in IC_FILES)


energy = load(ROOTS["Legacy"], "E_gam.txt")
spectra = {label: load_ic(root) for label, root in ROOTS.items()}
reference = spectra["Legacy"]

styles = {
    r"$2\pi$ only": ("tab:blue", "-"),
    "Lookup only": ("tab:orange", "-"),
    "Full corrected fork": ("tab:purple", "-"),
}

fig, axes = plt.subplots(3, 2, figsize=(9.4, 7.6), sharex=True, sharey=True)
for ax, (galaxy, index) in zip(axes.flat, GALAXIES):
    visible = reference[index] > np.nanmax(reference[index]) * 1.0e-8
    for label, (color, linestyle) in styles.items():
        ratio = np.divide(
            spectra[label][index],
            reference[index],
            out=np.full_like(reference[index], np.nan),
            where=visible,
        )
        ax.semilogx(energy, ratio, color=color, ls=linestyle, lw=1.7, label=label)
    ax.axhline(1.0, color="0.35", ls="--", lw=1.0)
    ax.axvline(ONE_KEV_GEV, color="0.5", ls=":", lw=1.0)
    ax.text(0.04, 0.92, galaxy, transform=ax.transAxes, va="top", fontweight="bold")
    ax.text(
        ONE_KEV_GEV, 0.97, "1 keV", transform=ax.get_xaxis_transform(),
        rotation=90, ha="right", va="top", color="0.35", fontsize=7,
    )
    ax.set_xlim(*X_LIMITS)
    ax.set_ylim(0.6, 2.0)
    ax.tick_params(which="both", direction="in", top=True, right=True)

for ax in axes[-1]:
    ax.set_xlabel(r"$E_\gamma\ [\mathrm{GeV}]$")
for ax in axes[:, 0]:
    ax.set_ylabel("IC / legacy IC")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
fig.suptitle("Independent effects of the two inverse-Compton corrections")
fig.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.12, hspace=0.08, wspace=0.08)
fig.savefig(DEST / "ic_isolated_corrections_ratio.png", dpi=300)
fig.savefig(DEST / "ic_isolated_corrections_ratio.pdf")


# A compact numerical summary at the energy where the reported dip occurs.
nearest = int(np.argmin(np.abs(np.log10(energy / ONE_KEV_GEV))))
rows = []
for galaxy, index in GALAXIES:
    values = [galaxy, energy[nearest]]
    values.extend(spectra[label][index, nearest] / reference[index, nearest] for label in styles)
    rows.append(values)

summary = DEST / "ic_isolated_corrections_1kev.csv"
np.savetxt(
    summary,
    np.asarray([row[1:] for row in rows]),
    delimiter=",",
    header="energy_GeV,ratio_2pi_only,ratio_lookup_only,ratio_full_corrected_fork",
    comments="",
    fmt="%.8e",
)
names = DEST / "ic_isolated_corrections_1kev_galaxies.txt"
names.write_text("\n".join(row[0] for row in rows) + "\n")

print(f"Saved {DEST / 'ic_isolated_corrections_ratio.png'}")
print(f"Saved {summary}")
for row in rows:
    print(f"{row[0]:9s}  2pi={row[2]:.4f}  lookup={row[3]:.4f}  both={row[4]:.4f}")
