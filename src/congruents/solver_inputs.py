"""Python-owned, pre-worker calculations from improvised legacy code/spectra.c.

No native callbacks: these arrays are complete before galaxy workers start.
Constants and algebra deliberately retain the scientific reference conventions.
"""
import math
import numpy as np
from scipy.special import hyp0f1
from . import constants as k
from .quadrature import integrate
from .serial import photon


def injection(energy, normalization, mass, cutoff):
    momentum = np.sqrt(energy**2 + 2*mass*energy)
    return normalization * momentum**(-2.2) * np.exp(-energy/cutoff) * momentum/energy


def normalization(mass, cutoff):
    # spectra_funcs.h:C_norm_E integrates kinetic energy from 0 to 1e8 GeV.
    def integrand(t):
        momentum = math.sqrt(t*t + 2*mass*t)
        return momentum**(-2.2) * math.exp(-t/cutoff) * momentum
    return integrate(integrand, 0., 1e8, rtol=1e-6)


def transport_inputs(rows, properties, kinetic):
    """Return transport[n,7,ne] and Cp[n]; secondary slot awaits native work.

    Preserve catalogue-index-10 suppression and the reference kinetic-energy
    convention for injection, even though the solver interpolates on total E.
    """
    n, ne = len(rows), len(kinetic)
    out = np.zeros((n, 7, ne), dtype=np.float64)
    norm_p = normalization(k.MP, 1e8)
    norm_e = normalization(k.ME, 1e5)
    cp = np.empty(n)
    for i, (row, p) in enumerate(zip(rows, properties)):
        h, nh, sigma, area, gas = p[0], p[1], p[3], p[4], p[5]
        power = row[3]*1.321680e-2*.1*1e51*k.ERG_GEV/k.YEAR
        ula = sigma/math.sqrt(2.)
        vai = 1000.*(ula/10.)/2.
        la = h/2.**3
        d0 = vai*la*1e5*k.PC
        collision = 1./(nh*1.4/1.17*40e-27*.5*k.C)
        loss = 1./(1./collision + 1./((h*k.PC)**2/d0))
        c = power*loss/(norm_p*2.*area*2.*h*k.PC**3)
        cp[i] = power*collision/norm_p
        def speed(mass, density, ion_fraction, cr_norm):
            momentum = np.sqrt(kinetic**2 + 2*mass*kinetic)**1.2
            return np.minimum(vai*(1.+2.3e-3*momentum*(density/1e3)**1.5*
                              ion_fraction*2./(ula/10.*cr_norm/2e-7)), k.C/1e5)
        vs = speed(k.MP, nh, 1., c)
        out[i,1] = vs*la*1e5*k.PC
        tau = 9.9*gas/1e3*h/1e2*1e27/out[i,1]
        gamma = 41.2*h/1e2*vs/1e3*1e27/out[i,1]
        out[i,0] = 1.-1./(hyp0f1(.25/1.25,tau/1.25**2) +
                            tau/gamma*hyp0f1(2.25/1.25,tau/1.25**2))
        if i == 10:
            out[i,0] *= .1
        out[i,2] = speed(k.ME, nh, 1., c)*la*1e5*k.PC
        out[i,3] = speed(k.ME, nh/1e3, 1./1e-4, (1.-out[i,0,0])*c)*la*1e5*k.PC
        out[i,4] = injection(kinetic, .2*power/norm_e, k.ME, 1e5)
        out[i,6] = injection(kinetic, cp[i], k.MP, 1e8)*out[i,0]
    if not np.all(np.isfinite(out)) or np.any(out < 0):
        raise RuntimeError("Invalid Python transport inputs")
    return out, cp


def free_free_inputs(rows, properties, energies):
    """Return [galaxy, 2, photon]: source FF (GeV^-1 s^-1), optical depth.

    Exact tau_FF_MK/eps_FF branches, fixed 1e4 K, including legacy subtraction
    and low-optical-depth series; no stabilization that changes the reference.
    """
    te = 1e4
    nu = energies/k.H
    gaunt = np.where(nu > 1e9,
        np.log(np.exp(5.960-math.sqrt(3)/math.pi*np.log(nu/1e9))+math.e),
        math.sqrt(3)/math.pi*(np.log((2*k.KB_ERG*te)**1.5 /
            (math.pi*k.E_ESU**2*math.sqrt(k.ME_G)*nu))-5*k.EULER/2))
    alpha = 4*k.E_ESU**6/(3*k.C*k.H_ERG)*math.sqrt(
        2*math.pi/(3*k.ME_G**3*k.KB_ERG*te)) * (
        1.-np.exp(-k.H_ERG*nu/(k.KB_ERG*te)))/nu**3*gaunt
    tau = 3*properties[:,6,None]/(k.YEAR*k.PC**2)*4.2e60*alpha/(8*2.54e-13)
    absorbed = np.where(tau > 1e-6, 1.-np.exp(-tau), tau-tau**2/2+tau**3/6)
    blackbody = np.array([photon("CMB", te, float(e)) for e in energies])
    ff = blackbody*k.C*4*math.pi*(rows[:,2,None]*1e3*k.PC)**2*absorbed
    result = np.ascontiguousarray(np.stack((ff,tau),axis=1))
    if not np.all(np.isfinite(result)) or np.any(result < 0):
        raise RuntimeError("Invalid Python free-free inputs")
    return result
