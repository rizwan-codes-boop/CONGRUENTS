"""Compare every standard legacy output file from two completed run directories.

Usage: python tests/compare_outputs.py C_OUTPUT PYTHON_OUTPUT
Exit nonzero on missing files, shape/NaN/zero/infinity mismatches or tolerance
failures. Reference-only debugging exports and Python-only provenance are ignored.
"""
import argparse
import json
from pathlib import Path
import numpy as np

FILES = ["T_CR.txt","E_gam.txt","fcal.txt","gal_data.txt","Urad_Ub.txt",
         "tau_ff.txt","E_loss_leptons.txt","E_loss_nucrit.txt","CR_specs.txt",
         "CR_specs_inj.txt","L_radio.txt","distmod.txt","L_gamma.txt",
         "spec_pi.txt","spec_pi_fcal1.txt","spec_nu.txt","spec_FF.txt"]
FILES += [f"spec_{p}_{population}_z{zone}.txt" for p in ("IC","BS","SY")
          for population in (1,2) for zone in (1,2) if p!="BS" or zone==1]
FILES += [f"tau_loss/tau_loss_z{zone}_{p}.txt" for zone in (1,2) for p in ("SY","BS","IC","DI","IO")]
FILES += [f"tau_loss/tau_loss_protons_{p}.txt" for p in ("PP","DI")]


def compare(reference, actual, tolerance=3e-6):
    report = {}
    for name in FILES:
        try:
            skip = 0 if name.startswith("CR_specs") else 1
            a = np.loadtxt(Path(reference)/name,skiprows=skip)
            b = np.loadtxt(Path(actual)/name,skiprows=skip)
            if a.shape!=b.shape or np.isnan(a).any() or np.isnan(b).any():
                raise ValueError("Shape mismatch or NaN")
            close = np.isclose(a,b,rtol=tolerance,atol=1e-280)
            # A reference zero must remain zero, even below the absolute floor.
            close &= (a!=0)|(b==0)
            finite = np.isfinite(a)&np.isfinite(b)&(a!=0)
            relative = np.abs((b[finite]-a[finite])/a[finite])
            report[name] = {"passed":bool(close.all()),"values":int(a.size),
                            "mismatches":int((~close).sum()),
                            "max_relative_error":float(relative.max()) if relative.size else 0.}
        except (OSError,ValueError) as exc:
            report[name] = {"passed":False,"error":str(exc)}
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("reference",type=Path)
    parser.add_argument("actual",type=Path)
    args = parser.parse_args()
    report = compare(args.reference,args.actual)
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if all(r["passed"] for r in report.values()) else 1)
