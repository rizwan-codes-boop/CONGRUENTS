"""Python loss-time, energy-budget and radio postprocessing of solver outputs."""
import bisect
import math
import numpy as np
from . import constants as k
from .serial import ionisation
from .quadrature import integrate
from .spectrum import integrate_spectrum
from .solver_inputs import free_free_inputs
from .preparation import PROPERTY_NAMES


def _lookup(table):
    x,y,z = table.x,table.y,table.values
    def value(a,b):
        if a<x[0] or a>x[-1] or b<y[0] or b>y[-1]:
            return 0.
        i,j = min(bisect.bisect_right(x,a)-1,len(x)-2),min(bisect.bisect_right(y,b)-1,len(y)-2)
        u,v = (a-x[i])/(x[i+1]-x[i]),(b-y[j])/(y[j+1]-y[j])
        return (1-v)*((1-u)*z[j,i]+u*z[j,i+1])+v*((1-u)*z[j+1,i]+u*z[j+1,i+1])
    return value


def _table_loss(energy, table_value, factor, tolerance):
    low,high = math.log(1e-3+k.ME),math.log(energy)
    if high<=low:
        return math.inf
    def integrand(x):
        final = math.exp(x)
        return final*(energy-final)*factor*table_value(energy-final,energy) if final<energy else 0.
    rate = integrate(integrand,low,high,rtol=tolerance)
    return energy/rate if rate else math.inf


def _sync(energy,b):
    gamma = energy/k.ME
    return energy/(k.SIGMA_MB*k.MB_CM2*k.C*gamma**2*(1-1/gamma**2)*b**2/(8*math.pi)*k.ERG_GEV)


def _plasma(energy,nh):
    plasma = math.sqrt(nh*1.1/(math.pi*k.ME_G))*k.E_ESU
    rate = .75*k.C*k.SIGMA_MB*k.MB_CM2*k.ME*nh*1.1*(math.log(energy/k.ME)+2*math.log(k.ME/(k.H*plasma)))
    return energy/rate


def calculate(preparation,kinetic,electron,photon,transport,electrons,emission,radio,gamma,bs):
    n,ne = len(transport),len(electron)
    out = {}
    rows = np.array([g.as_row() for g in preparation.galaxies])
    props = np.array([[p[key] for key in PROPERTY_NAMES] for p in preparation.properties])
    out["gal_data"] = props[:,:9].copy()
    out["Urad_Ub"] = np.array([preparation.radiation(i,[])['urad_ub_ev_cm3'] for i in range(n)])
    escape,critical = np.empty((n,2,ne)),np.empty((n,10))
    for zone in (1,2):
        for process in ("SY","BS","IC","DI","IO"):
            out[f"tau_loss_z{zone}_{process}"] = np.empty((n,ne))
    out["tau_loss_protons_PP"] = np.empty((n,ne))
    out["tau_loss_protons_DI"] = np.empty((n,ne))
    brems = _lookup(bs)
    for i,p in enumerate(props):
        ic = _lookup(gamma[i])
        for zone,nh,b,h,d in ((1,p[1],p[2],p[0],transport[i,2]),
                               (2,p[1]/1000,p[9],50*p[0],transport[i,3])):
            out[f"tau_loss_z{zone}_SY"][i] = [_sync(e,b) for e in electron]
            out[f"tau_loss_z{zone}_BS"][i] = [_table_loss(e,brems,k.C*nh*k.MB_CM2,1e-6) for e in electron]
            out[f"tau_loss_z{zone}_IC"][i] = [_table_loss(e,ic,1.,1e-8) for e in electron]
            out[f"tau_loss_z{zone}_DI"][i] = (h*k.PC)**2/d
            # Reference time-grid diagnostic uses ionisation in both zones;
            # its halo solver and critical-energy diagnostic use plasma losses.
            out[f"tau_loss_z{zone}_IO"][i] = -electron/np.array(ionisation(electron,nh))
            ec = math.sqrt(2*1.49e9*k.ME_G*k.C/(3*b*k.E_ESU))*math.pi*k.ME
            values = (_table_loss(ec,brems,k.C*nh*k.MB_CM2,1e-6),_sync(ec,b),
                      _table_loss(ec,ic,1.,1e-8),
                      -ec/ionisation([ec],nh)[0] if zone==1 else _plasma(ec,nh),
                      (h*k.PC)**2/np.interp(ec,electron,d))
            critical[i,5*(zone-1):5*zone] = ec/np.array(values)
        escape[i] = electrons[i,:2]/out["tau_loss_z1_DI"][i]
        out["tau_loss_protons_PP"][i] = 1/(p[1]*1.4/1.17*40e-27*.5*k.C)
        out["tau_loss_protons_DI"][i] = (p[0]*k.PC)**2/transport[i,1]
    out["E_loss_nucrit"] = critical
    budget = np.zeros((n,16))
    for i in range(n):
        budget[i,:2] = [integrate_spectrum(kinetic,transport[i,j]) for j in (4,5)]
        for population,indices,start in ((0,(4,0,2,8,6),2),(1,(5,1,3,9,7),9)):
            budget[i,start:start+5] = [integrate_spectrum(photon,emission[i,j])/budget[i,population] for j in indices]
            budget[i,start+6] = integrate_spectrum(kinetic,escape[i,population])/budget[i,population]
    out["E_loss_leptons"] = budget
    out["CR_specs"] = np.concatenate((transport[:,6:7],electrons),axis=1).reshape(n*5,ne)
    out["CR_specs_inj"] = np.concatenate((transport[:,4:6],escape),axis=1).reshape(n*4,ne)
    ff = free_free_inputs(rows,props,np.array([1.49e9*k.H]))
    radio = radio.copy()
    radio[:,:2] *= np.exp(-ff[:,1,:])
    out["L_radio"] = np.column_stack((radio,ff[:,0,0]))*(1.49e9*6.62607015e-34*k.H)
    for name,a in out.items():
        # Infinite time is meaningful where a tabulated loss rate is zero.
        if np.isnan(a).any() or np.any(a<0) or (not name.startswith("tau_loss") and not np.isfinite(a).all()):
            raise RuntimeError(f"Invalid diagnostic: {name}")
        a.flags.writeable = False
    return out
