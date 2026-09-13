"""Serial Python attenuation, retaining improvised legacy code table and kernel rules."""
import hashlib
import math
from pathlib import Path
import numpy as np
from . import constants as k
from .serial import dilution, photon
from .quadrature import integrate


class EBLTable:
    def __init__(self, path):
        self.path = Path(path)
        raw = self.path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        tokens = [t for line in raw.decode().splitlines()
                  if not line.lstrip().startswith("#") for t in line.split()]
        if len(tokens)<2:
            raise ValueError("Missing EBL dimensions")
        nz, ne = int(tokens[0]), int(tokens[1])
        if nz < 2 or ne < 2 or len(tokens) != 2+nz+ne+nz*ne:
            raise ValueError("Invalid EBL table dimensions/token count")
        data = np.array(tokens[2:], dtype=float)
        self.z, self.energy = data[:nz], data[nz:nz+ne]
        self.values = data[nz+ne:].reshape(nz,ne)
        if (not np.isfinite(data).all() or np.any(np.diff(self.z)<=0) or
                np.any(np.diff(self.energy)<=0) or np.any(self.energy<=0) or
                np.any(self.z<0) or np.any(self.values<0)):
            raise ValueError("Invalid EBL coordinates/values")
        for a in (self.z,self.energy,self.values):
            a.flags.writeable = False

    def optical_depth(self, source_energy_gev, redshifts):
        """Bilinear extrapolation, E/(1+z) in eV capped at 1e15; clip tau>=0.

        This deliberately preserves spectra.c's source-grid indexing rather
        than resampling optical depths alongside the emitted spectrum.
        """
        e, z = np.asarray(source_energy_gev,float), np.asarray(redshifts,float)
        if (e.ndim!=1 or z.ndim!=1 or not np.isfinite(e).all() or
                not np.isfinite(z).all() or np.any(e<=0) or np.any(z<0)):
            raise ValueError("Invalid EBL query")
        x = np.minimum(e[None,:]/(1+z[:,None])*1e9,1e15)
        ix = np.clip(np.searchsorted(self.energy,x,side="right")-1,0,len(self.energy)-2)
        iy = np.clip(np.searchsorted(self.z,z,side="right")-1,0,len(self.z)-2)[:,None]
        u = (x-self.energy[ix])/(self.energy[ix+1]-self.energy[ix])
        v = (z[:,None]-self.z[iy])/(self.z[iy+1]-self.z[iy])
        a = self.values
        return np.maximum(0.,(1-u)*(1-v)*a[iy,ix]+u*(1-v)*a[iy,ix+1]+
                          (1-u)*v*a[iy+1,ix]+u*v*a[iy+1,ix+1])


def internal_optical_depth(preparation, energies):
    """tau_gg_gal_BW: original log-target-energy quadrature and path length."""
    e = np.asarray(energies,float)
    if e.ndim!=1 or not np.isfinite(e).all() or np.any(e<=0):
        raise ValueError("Invalid gamma-ray energy grid")
    out = np.zeros((len(preparation.galaxies),len(e)))
    with preparation._lock:
        preparation._ensure_open()
        low, high = map(math.log,preparation.config[4:6])
        for i,(galaxy,p) in enumerate(zip(preparation.galaxies,preparation.properties)):
            dust = p["dust_temperature_k"]
            weights = dilution(galaxy,dust)
            def density(target):
                return (photon("CMB",k.TCMB*(1+galaxy.redshift),target)+
                        sum(w*photon(f,dust,target) for w,f in
                            zip(weights,("3000","4000","7500","UV","FIR"))))
            for j,energy in enumerate(e):
                def integrand(x):
                    target = math.exp(x)
                    product = energy*target
                    if product < k.ME**2:
                        return 0.
                    beta = math.sqrt(1-k.ME**2/product)
                    cross = 3/16*k.SIGMA_MB*(1-beta**2)*(2*beta*(beta**2-2)+
                             (3-beta**4)*math.log((1+beta)/(1-beta)))
                    return target*density(target)*cross
                out[i,j] = p["height_pc"]*k.PC*k.MB_CM2*integrate(integrand,low,high,rtol=1e-6)
    return out
