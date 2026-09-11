"""Print a reduced-grid reference driver; do not import the Python/C adapter.

Only resolution constants and extra raw-output statements change. All physics,
loop bodies and solver calls come from the preserved production source.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
kind = sys.argv[1]
source = (ROOT/("spectra.c" if kind == "spectra" else "create_interp_objects.c")).read_text()


def replace_once(old, new, count=1):
    global source
    if source.count(old) != count:
        raise RuntimeError(f"Reference source changed: expected {count} occurrences of {old!r}")
    source = source.replace(old, new)


medium = "--medium" in sys.argv
cells, photons = (32,16) if medium else (16,8)
replace_once("size_t n_pts[2] = { 1000, 1000 };",
             f"size_t n_pts[2] = {{ {cells}, {cells} }};")
if kind == "spectra":
    replace_once("unsigned int n_T_CR = 1000;", f"unsigned int n_T_CR = {cells};")
    replace_once("int n_E_gam = 500;", f"int n_E_gam = {photons};")
    replace_once("E_CRe_lims__GeV, 500,", f"E_CRe_lims__GeV, {cells},", 2)
    exports = {
        "fcal": ("f_cal", "n_T_CR"),
        "Dp": ("D__cm2sm1", "n_T_CR"),
        "Dd": ("D_e__cm2sm1", "n_T_CR"),
        "Dh": ("D_e_z2__cm2sm1", "n_T_CR"),
        "Q1": ("Q_e_1_z1", "n_T_CR"),
        "Q2": ("Q_e_2_z1", "n_T_CR"),
        "protons": ("q_p_SS_z1", "n_T_CR"),
    }
    for name, array in zip(("e1d","e2d","e1h","e2h"),
                          ("q_e_SS_1_z1","q_e_SS_2_z1","q_e_SS_1_z2","q_e_SS_2_z2")):
        exports[name] = (array, "n_T_CR")
    for name, array in zip(("ic1d","ic2d","bs1d","bs2d","sy1d","sy2d","ic1h","ic2h",
                            "sy1h","sy2h","ff","tau_ff","pi","pi_fcal1","nu"),
                          ("spec_IC_1_z1","spec_IC_2_z1","spec_BS_1_z1","spec_BS_2_z1",
                           "spec_SY_1_z1","spec_SY_2_z1","spec_IC_1_z2","spec_IC_2_z2",
                           "spec_SY_1_z2","spec_SY_2_z2","spec_FF","tau_FF","spec_pi",
                           "spec_pi_fcal1","spec_nu")):
        exports[name] = (array, "n_E_gam")
    writes = "\n".join(f'write_2D_file(n_gal,{size},{array},"reference",'
                       f'string_cat(outfp,"/{name}.txt"));'
                       for name,(array,size) in exports.items())
    anchor = "    double **array2Dlist[4]"
    # Stop at the Week-3 boundary; observer/output postprocessing is Week 4.
    replace_once(anchor, writes+"\nreturn 0;\n"+anchor)
if kind == "precompute" and "--full-precision" in sys.argv:
    # Match Python's full-precision in-memory tables, not the legacy text
    # writer's six-digit rounding. Otherwise its truncated upper bounds can
    # exclude an endpoint and turn a nonzero interpolated value into zero.
    prefix = r'''
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
static int reference_fprintf(FILE *stream,const char *format,...) {
    char expanded[4096]; size_t j=0;
    for(size_t i=0;format[i] && j<4080;i++) {
        if(!strncmp(format+i,"%le",3)) {
            memcpy(expanded+j,"%.17e",5);j+=5;i+=2;
        } else if(!strncmp(format+i,"%e",2)) {
            memcpy(expanded+j,"%.17e",5);j+=5;i++;
        } else expanded[j++]=format[i];
    }
    expanded[j]=0;
    va_list args;va_start(args,format);
    int result=vfprintf(stream,expanded,args);
    va_end(args);return result;
}
#define fprintf reference_fprintf
'''
    source = prefix + source
print(source)
