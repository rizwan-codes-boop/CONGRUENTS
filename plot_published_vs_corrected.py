"""Compare archived paper-production spectra with corrected CONGRUENTS-i outputs.

Produces a six-panel IC overlay and a six-panel component-ratio figure.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
NEW = HERE / "output"
PUBLISHED = HERE.parents[1] / "CONGRuENTS_I_files" / "local_final"

GALAXIES = [
    ("NGC 2146", 4),
    ("NGC 2403", 3),
    ("M82", 5),
    ("NGC 253", 2),
    ("SMC", 9),
    ("M33", 7),
]

PLOT_LIMITS = (1.0e-16, 1.0e5)
ONE_KEV_GEV = 1.0e-6


def load(root, filename):
    return np.loadtxt(root / filename, skiprows=1)


def combined(root, filenames):
    return sum(load(root, filename) for filename in filenames)


energy = load(NEW, "E_gam.txt")

FILES = {
    "IC": [
        "spec_IC_1_z1.txt", "spec_IC_2_z1.txt",
        "spec_IC_1_z2.txt", "spec_IC_2_z2.txt",
    ],
    "Synchrotron": [
        "spec_SY_1_z1.txt", "spec_SY_2_z1.txt",
        "spec_SY_1_z2.txt", "spec_SY_2_z2.txt",
    ],
    "Bremsstrahlung": ["spec_BS_1_z1.txt", "spec_BS_2_z1.txt"],
    "Free-free": ["spec_FF.txt"],
    "Pion": ["spec_pi.txt"],
}

published = {name: combined(PUBLISHED, files) for name, files in FILES.items()}
corrected = {name: combined(NEW, files) for name, files in FILES.items()}
published["Total"] = sum(published.values())
corrected["Total"] = sum(corrected.values())


def significant_ratio(new, reference, fraction=1.0e-8):
    """Return a ratio only where the reference curve is scientifically visible."""
    threshold = np.nanmax(reference) * fraction
    return np.divide(
        new,
        reference,
        out=np.full_like(new, np.nan),
        where=reference > threshold,
    )


def style_axis(ax, name):
    ax.set_xscale("log")
    ax.set_xlim(*PLOT_LIMITS)
    ax.axvline(ONE_KEV_GEV, color="0.45", ls=":", lw=1)
    ax.text(
        ONE_KEV_GEV, 0.97, "1 keV", transform=ax.get_xaxis_transform(),
        rotation=90, ha="right", va="top", color="0.35", fontsize=7,
    )
    ax.text(0.04, 0.92, name, transform=ax.transAxes, va="top", fontweight="bold")
    ax.tick_params(which="both", direction="in", top=True, right=True)


# Figure 1: direct IC comparison. Published is dashed; corrected is solid.
fig, axes = plt.subplots(3, 2, figsize=(9.2, 7.5), sharex=True, sharey=False)
for ax, (name, index) in zip(axes.flat, GALAXIES):
    ax.loglog(
        energy, published["IC"][index], color="tab:purple", ls="--", lw=1.7,
        label="Published-production IC",
    )
    ax.loglog(
        energy, corrected["IC"][index], color="tab:purple", ls="-", lw=2.0,
        label="Corrected IC",
    )
    visible_x = (energy >= PLOT_LIMITS[0]) & (energy <= PLOT_LIMITS[1])
    pair = np.concatenate((published["IC"][index, visible_x], corrected["IC"][index, visible_x]))
    peak = np.nanmax(pair)
    significant = pair[pair > peak * 1.0e-8]
    ax.set_ylim(np.nanmin(significant) * 0.5, peak * 2.0)
    style_axis(ax, name)

for ax in axes[-1]:
    ax.set_xlabel(r"$E_\gamma\ [\mathrm{GeV}]$")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$E_\gamma^2\phi_\gamma\ [\mathrm{GeV\,s^{-1}\,cm^{-2}}]$")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)
fig.suptitle("Inverse-Compton spectrum: published production vs corrected model")
fig.subplots_adjust(left=0.11, right=0.98, top=0.93, bottom=0.11, hspace=0.08, wspace=0.08)
fig.savefig(NEW / "published_vs_corrected_ic.png", dpi=300)
fig.savefig(NEW / "published_vs_corrected_ic.pdf")


# Figure 2: corrected / published ratios for all major components.
colors = {
    "Total": "black",
    "IC": "tab:purple",
    "Synchrotron": "tab:brown",
    "Bremsstrahlung": "#b5a800",
    "Free-free": "orangered",
    "Pion": "tab:blue",
}
fig, axes = plt.subplots(3, 2, figsize=(9.2, 7.5), sharex=True, sharey=True)
for ax, (name, index) in zip(axes.flat, GALAXIES):
    for component in ("Total", "IC", "Synchrotron", "Bremsstrahlung", "Free-free", "Pion"):
        ratio = significant_ratio(corrected[component][index], published[component][index])
        ax.semilogx(energy, ratio, color=colors[component], lw=1.5, label=component)
    ax.axhline(1.0, color="0.35", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_ylim(6.0e-1, 2.0e0)
    style_axis(ax, name)

for ax in axes[-1]:
    ax.set_xlabel(r"$E_\gamma\ [\mathrm{GeV}]$")
for ax in axes[:, 0]:
    ax.set_ylabel("Corrected / published")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
fig.suptitle("Spectral changes after implementation corrections")
fig.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.13, hspace=0.08, wspace=0.08)
fig.savefig(NEW / "published_vs_corrected_ratios.png", dpi=300)
fig.savefig(NEW / "published_vs_corrected_ratios.pdf")

plt.show()
