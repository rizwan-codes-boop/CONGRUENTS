"""Serial model preparation in Python, preserving improvised legacy code equations."""
import math
import sys
import numpy as np
from scipy.special import kv
from . import constants as k
from .quadrature import integrate

def log_grid(n, low, high):
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("Grid size must be an integer >= 2")
    if not 0 < low < high or not math.isfinite(high):
        raise ValueError("Grid bounds must be finite, positive and increasing")
    lo = math.log(low)
    step = (math.log(high)-lo)/(n-1)
    return np.array([math.exp(lo+step*i) for i in range(n)], dtype=np.float64)

def ionisation(energies, density):
    if not math.isfinite(density) or density < 0:
        raise ValueError("density must be finite and non-negative")
    out = []
    for e in energies:
        if not math.isfinite(e) or e <= k.ME:
            raise ValueError("Expected finite TOTAL electron energy above rest mass")
        v = -9./4.*k.C*k.SIGMA_MB*k.MB_CM2*k.ME*density*1.1*(
            math.log(e/k.ME)+.91*2./3.*math.log(k.ME*1e9/15.)+
            2.*.09*2./3.*math.log(k.ME*1e9/41.5))
        if not math.isfinite(v):
            raise RuntimeError("Non-finite ionisation result")
        out.append(v)
    return out

def bounds(galaxies, properties):
    maxb = max(p["magnetic_field_gauss"] for p in properties)
    minb = min(p["halo_magnetic_field_gauss"] for p in properties)
    maxz = max(g.redshift for g in galaxies)
    lowe, highe = 1e-3+k.ME, 1e8+k.ME
    common = 2*math.pi**2*k.ME_G**2*k.C**3
    sx0 = common/(3*k.E_ESU*maxb*k.H_ERG)*(1e-16*k.ME)/highe**2*.1
    sx1 = common/(3*k.E_ESU*minb/10*k.H_ERG)*(1e8*(1+maxz)*k.ME)/lowe**2*10
    config = (1e-16, 1e8, lowe, highe, 1.59362*k.KB*k.TCMB*1e-4, 1e-7, sx0, sx1)
    ts = (k.TCMB*(1+maxz), min(p["dust_temperature_k"] for p in properties),
          max(p["dust_temperature_k"] for p in properties))
    return config, ts

def temperature_grid(field, temperatures):
    if field == "CMB":
        low, high, step = k.TCMB, temperatures[0], .5
        limits = (math.floor(k.TCMB*10)/10, math.ceil(high*10)/10)
    elif field == "FIR":
        low, high, step = temperatures[1], temperatures[2], 5.
        limits = (math.floor(low), math.ceil(high))
    else:
        raise ValueError("Temperature grid is only for CMB/FIR")
    n = math.ceil((high-low)/step)+1
    if n < 2 or n > 4096 or limits[0] <= 0 or limits[1] <= limits[0]:
        raise ValueError("Legacy temperature grid requires at least two planes")
    delta = (limits[1]-limits[0])/(n-1)
    return tuple([limits[0]] + [limits[0]+i*delta for i in range(1,n-1)] + [limits[1]])

def photon(field, temp, energy):
    if field == "UV":
        wavelength = k.C*k.H/energy*1e4
        if .134 < wavelength <= .246:
            return 2.373/energy**2*k.ERG_GEV*wavelength**(-.6678)
        if .110 < wavelength <= .134:
            return 68.25/energy**2*k.ERG_GEV*wavelength
        if .0912 < wavelength <= .110:
            return 1.287e5/energy**2*k.ERG_GEV*wavelength**4.4172
        return 0.
    if field in ("3000","4000","7500"):
        temp = float(field)
    exponent = energy/(k.KB*temp)
    if exponent > math.log(sys.float_info.max):
        return 0.
    value = 8*math.pi*energy**2/(k.H*k.C)**3/(math.exp(exponent)-1)
    return value*energy/(2e12*k.H) if field == "FIR" else value

def dilution(galaxy, dust):
    mass, radius, sfr = galaxy.stellar_mass_msun, galaxy.radius_kpc, galaxy.sfr_msun_per_year
    old = 10**(.8480565*math.log10(mass/.56)+1.521623)
    young = (10**(.7969616*math.log10(sfr/.56)+9.007323) if math.log10(sfr)>-2.6
             else 10**(.9867204*math.log10(sfr/.56)+9.476592))
    fir = 10**(1.096548*math.log10(sfr/.56)+9.710084)
    densities = (k.ARAD*3000**4,k.ARAD*4000**4,k.ARAD*7500**4,4450.1668,
                 24.8863*8*math.pi/((k.H*k.C)**3*2e12*k.H)*(k.KB*dust)**5)
    luminosities = (.574*old,.426*old,.763*young,.237*young,fir)
    return tuple(l*k.LSUN/(u*2*math.pi*(radius*1e3*k.PC)**2*k.C)
                 for l,u in zip(luminosities,densities))

def radiation(galaxy, properties, energies, config):
    dust, cmb = properties["dust_temperature_k"], k.TCMB*(1+galaxy.redshift)
    d = dilution(galaxy,dust)
    def components(e):
        c = photon("CMB",cmb,e)
        f = d[4]*photon("FIR",dust,e)
        a,b,s,v = (d[i]*photon(field,0,e) for i,field in enumerate(("3000","4000","7500","UV")))
        # Legacy diagnostic helpers use unsplit old luminosity.
        return (c,f,a/.574,b/.426,s,v,c+f+a+b+s+v)
    energies = tuple(float(x) for x in energies)
    if any(not math.isfinite(e) or e<=0 for e in energies):
        raise ValueError("Photon energies must be finite and positive")
    fields = tuple(zip(*(components(e) for e in energies))) if energies else ((),)*7
    integrals = [integrate(lambda x: math.exp(x)**2*components(math.exp(x))[j],
                           math.log(config[4]),math.log(config[5]),rtol=1e-6)*1e9 for j in range(7)]
    magnetic = properties["magnetic_field_gauss"]**2/(8*math.pi)/k.GEV_ERG*1e9
    return fields, (integrals[6],magnetic,*integrals[:6])

def _ic_kernel(electron, target, emitted, density):
    if electron == emitted:
        return 0.
    gamma = 4*target/k.ME*electron/k.ME
    q = emitted/(gamma*(electron-emitted))
    if (k.ME/(2*electron))**2 < q < 1:
        return density*(2*q*math.log(q)+(1+2*q)*(1-q)+(gamma*q)**2*(1-q)/(2*(1+gamma*q)))
    return 0.

def generate(kind, field, temperature, nx, ny, config):
    from .inputs import Grid
    Grid(nx, ny)
    if kind not in ("emission","gamma","bs","sy"):
        raise ValueError("Unknown table kind")
    if field not in ("3000", "4000", "7500", "UV", "CMB", "FIR"):
        raise ValueError("Unknown radiation field")
    if field in ("CMB","FIR") and (not math.isfinite(temperature) or not 2<=temperature<=1000):
        raise ValueError("Invalid field temperature")
    if len(config)!=8 or any(not math.isfinite(v) or v<=0 for v in config):
        raise ValueError("Invalid table bounds")
    if any(config[i]>=config[i+1] for i in (0,2,4,6)):
        raise ValueError("Table bounds must increase")
    if kind == "sy":
        x = log_grid(nx,config[6],config[7])
        end = min(math.log(config[7]),math.log(150.))
        values = [a*integrate(lambda z: math.exp(z)*float(kv(5./3.,math.exp(z))),
                             math.log(a),end) if math.log(a)<end else 0. for a in x]
        return x,np.array([1.]),np.array(values).reshape(1,nx)
    x = log_grid(nx,config[4] if kind=="gamma" else config[0],
                 config[3] if kind=="gamma" else config[1])
    y = log_grid(ny,config[2],config[3])
    z = np.zeros((ny,nx))
    delta_axis = [0,.01,.02,.05,.1,.2,.5,1,2,5,10]
    phi1 = [45.79,45.43,45.09,44.11,42.64,40.16,34.97,29.97,24.73,18.09,13.65]
    phi2 = [44.46,44.38,44.24,43.65,42.49,40.19,34.93,29.78,24.34,17.28,12.41]
    for i, out in enumerate(x):
        for j, electron in enumerate(y):
            out,electron = float(out),float(electron)
            if kind == "bs":
                if out >= electron:
                    continue
                delta = out*k.ME/(4*k.ALPHA*electron*(electron-out))
                if delta > 2:
                    a=b=4*(math.log(2*electron/k.ME*((electron-out)/out))-.5)
                else:
                    a,b = float(np.interp(delta,delta_axis,phi1)),float(np.interp(delta,delta_axis,phi2))
                z[j,i] = 3/(8*math.pi)*k.SIGMA_MB*k.ALPHA*max(
                    (1+(1-out/electron)**2)*a-2/3*(1-out/electron)*b,0)/out
                continue
            if kind == "emission":
                # C fmax/fmin ignore NaN bounds when photon >= electron;
                # the kernel itself is zero there.
                if out >= electron:
                    continue
                low=max(math.log(config[4]),math.log(out*k.ME**2/(4*electron*(electron-out))))
                high=min(math.log(config[5]),math.log(out*electron/(electron-out)))
            else:
                if out >= electron:
                    continue
                low,high = math.log(config[4]),math.log(config[5])
            if low < high:
                def integrand(logtarget):
                    target=math.exp(logtarget)
                    emitted=out+target if kind=="gamma" else out
                    return _ic_kernel(electron,target,emitted,photon(field,temperature,target))
                result=integrate(integrand,low,high)
                z[j,i]=.75*k.SIGMA_MB*k.MB_CM2*k.C/(electron/k.ME)**2*max(result,0)
    return x,y,z
