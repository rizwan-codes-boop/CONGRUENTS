"""Week-3 galaxy solver: Python orchestration, native galaxy-loop kernels."""
import ctypes as ct
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
import sys
import threading
import numpy as np
from . import constants, serial
from .solver_inputs import transport_inputs, free_free_inputs
from .preparation import Preparation, PROPERTY_NAMES, Table, TableData

_LOCK = threading.RLock()
_DP = ct.POINTER(ct.c_double)
_IP = ct.POINTER(ct.c_int)


class _TableInput(ct.Structure):
    _fields_ = [("nx", ct.c_size_t), ("ny", ct.c_size_t),
                ("x", _DP), ("y", _DP), ("values", _DP)]


def _pointer(a):
    return a.ctypes.data_as(_DP)


def _descriptor(table):
    # Preparation retains the arrays until the synchronous batch has returned.
    return _TableInput(len(table.x), len(table.y), _pointer(table.x),
                       _pointer(table.y), _pointer(table.values))


def _load(path=None):
    suffix = "dylib" if sys.platform == "darwin" else "so"
    path = Path(path) if path else Path(__file__).resolve().parents[2]/"build"/f"libcongruents_solver.{suffix}"
    if not path.is_file():
        raise FileNotFoundError(f"Build the optional solver first: make solver ({path})")
    lib = ct.CDLL(str(path.resolve()))
    lib.cg_solver_abi.argtypes = []
    lib.cg_solver_abi.restype = ct.c_uint
    if lib.cg_solver_abi() != 2:
        raise RuntimeError("Unsupported solver ABI")
    lib.cg_transport_batch.argtypes = [ct.c_int, ct.c_size_t, ct.c_size_t,
                                      _DP, _DP, _DP, _DP, _DP, _IP]
    lib.cg_transport_batch.restype = ct.c_int
    lib.cg_solver_batch.argtypes = [ct.c_int, *([ct.c_size_t]*4),
                                   *([_DP]*7), *([ct.POINTER(_TableInput)]*4),
                                   _DP, _DP, _IP]
    lib.cg_solver_batch.restype = ct.c_int
    return lib


def _check(code, statuses):
    if code:
        raise RuntimeError(f"Native batch rejected input or resources (status {code})")
    failed = np.flatnonzero(statuses)
    if failed.size:
        raise RuntimeError(f"Native calculation failed for catalogue indices {failed.tolist()}: "
                           f"statuses {statuses[failed].tolist()}; outputs discarded")


def _readonly(a):
    a.flags.writeable = False
    return a


@dataclass(frozen=True)
class SolverGrid:
    cosmic_rays: int = 1000
    photons: int = 500
    cells: int = 500

    def __post_init__(self):
        for n, maximum in ((self.cosmic_rays, 4096), (self.photons, 4096), (self.cells, 500)):
            if isinstance(n, bool) or not isinstance(n, int):
                raise TypeError("Solver grid sizes must be integers")
            if not 4 <= n <= maximum:
                raise ValueError(f"Solver grid size must be in [4, {maximum}]")


@dataclass(frozen=True)
class SolverResult:
    kinetic_energy_gev: np.ndarray
    electron_energy_gev: np.ndarray
    photon_energy_gev: np.ndarray
    transport: dict
    electrons: dict
    emission: dict
    metadata: dict

    def save(self, path):
        """Save named arrays without pickle; this is not legacy output format."""
        arrays = {"kinetic_energy_gev": self.kinetic_energy_gev,
                  "electron_energy_gev": self.electron_energy_gev,
                  "photon_energy_gev": self.photon_energy_gev,
                  "metadata_json": np.array(json.dumps(self.metadata, sort_keys=True))}
        for group in ("transport", "electrons", "emission"):
            arrays.update({f"{group}__{name}": value for name, value in getattr(self, group).items()})
        np.savez_compressed(path, **arrays)


TRANSPORT = ("fcal", "proton_diffusion_cm2_s", "disc_diffusion_cm2_s",
             "halo_diffusion_cm2_s", "primary_injection_gev_s", "secondary_injection_gev_s",
             "proton_steady_state_gev")
ELECTRONS = ("primary_disc", "secondary_disc", "primary_halo", "secondary_halo")
EMISSION = ("IC_primary_disc", "IC_secondary_disc", "BS_primary_disc", "BS_secondary_disc",
            "SY_primary_disc", "SY_secondary_disc", "IC_primary_halo", "IC_secondary_halo",
            "SY_primary_halo", "SY_secondary_halo", "free_free", "tau_free_free",
            "pion", "pion_full_calorimetry", "neutrino")


def solve(preparation, grid=None, threads=1, library=None, *, legacy_table_precision=True):
    """Return source components for all catalogue rows, in order.

    Electron outputs: GeV^-1; emission: GeV^-1 s^-1, except dimensionless
    tau_free_free. Disc components include baseline free-free attenuation.
    This does not yet apply EBL or gamma-gamma attenuation or observer distances.
    Default table rounding reproduces the legacy %.6e cache-file round trip.
    Set legacy_table_precision=False only for a full-precision comparison.
    """
    if not isinstance(preparation, Preparation):
        raise TypeError("Expected Preparation")
    grid = SolverGrid() if grid is None else grid
    if not isinstance(grid, SolverGrid):
        raise TypeError("Expected SolverGrid")
    if not isinstance(legacy_table_precision, bool):
        raise TypeError("legacy_table_precision must be bool")
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 2147483647:
        raise ValueError("threads must be a positive C integer")
    with preparation._lock, _LOCK:
        preparation._ensure_open()
        lib = _load(library)
        n, ne, np_ = len(preparation.galaxies), grid.cosmic_rays, grid.photons
        kinetic = serial.log_grid(ne, 1e-3, 1e8)
        electron = serial.log_grid(ne, 1e-3+constants.ME, 1e8+constants.ME)
        photon = serial.log_grid(np_, 1e-16, 1e8)
        rows = np.array([g.as_row() for g in preparation.galaxies], dtype=np.float64)
        properties = np.array([[p[k] for k in PROPERTY_NAMES] for p in preparation.properties])
        transport, cp = transport_inputs(rows, properties, kinetic)
        free_free = free_free_inputs(rows, properties, photon)
        density = np.ascontiguousarray(properties[:,1])
        fcal = np.ascontiguousarray(transport[:,0,:])
        secondary = np.zeros((n,ne))
        statuses = np.zeros(n, dtype=np.int32)
        code = lib.cg_transport_batch(threads, n, ne, _pointer(kinetic), _pointer(density),
                _pointer(cp), _pointer(fcal), _pointer(secondary), statuses.ctypes.data_as(_IP))
        _check(code, statuses)
        transport[:,5,:] = secondary
        props = np.ascontiguousarray(np.column_stack(
            (properties[:,0], properties[:,1], properties[:,2], properties[:,9],
             properties[:,6], rows[:,2], rows[:,3], cp)))
        diffusion = np.ascontiguousarray(transport[:,2:4,:])
        injection = np.ascontiguousarray(transport[:,4:6,:])
        electrons = np.zeros((n,4,ne))
        emission = np.zeros((n,15,np_))
        fcal = np.ascontiguousarray(transport[:,0,:])
        combined = []
        rounded = {}
        def round_values(values):
            a = np.asarray(values)
            return np.array([float(format(float(v), ".6e")) for v in a.ravel()]).reshape(a.shape)
        def temperatures(field):
            values = preparation.temperature_grid(field)
            return tuple(round_values(values)) if legacy_table_precision else values
        def table_at(kind, field="3000", temperature=0.):
            key = (kind, field, temperature)
            # The legacy loader checks Gamma bounds against emission bounds,
            # rejects those files and regenerates Gamma at full precision.
            if not legacy_table_precision or kind == "gamma":
                return preparation.table(kind,field,temperature)
            if key not in rounded:
                actual = temperature
                if field in ("CMB","FIR"):
                    nodes = temperatures(field)
                    actual = preparation.temperature_grid(field)[nodes.index(temperature)]
                table = preparation.table(kind,field,actual)
                rounded[key] = Table(TableData(round_values(table.x), round_values(table.y),
                                                round_values(table.values)), table.metadata)
            return rounded[key]
        try:
            for kind in ("emission", "gamma"):
                family = []
                combined.append(family)
                for i in range(n):
                    family.append(preparation.combined_ic(i, kind, _table_provider=table_at,
                        _temperature_provider=(preparation.temperature_grid if kind == "gamma" else temperatures)))
            ic = (_TableInput*n)(*map(_descriptor, combined[0]))
            gamma = (_TableInput*n)(*map(_descriptor, combined[1]))
            bs = _descriptor(table_at("bs"))
            sy = _descriptor(table_at("sy"))
            code = lib.cg_solver_batch(threads,n,ne,np_,grid.cells,
                _pointer(electron),_pointer(photon),_pointer(props),_pointer(diffusion),_pointer(injection),
                _pointer(kinetic),_pointer(fcal),
                ic,gamma,ct.byref(bs),ct.byref(sy),_pointer(electrons),_pointer(emission),
                statuses.ctypes.data_as(_IP))
            _check(code,statuses)
            emission[:,:6,:] *= np.exp(-free_free[:,1,None,:])
            emission[:,10:12,:] = free_free
        finally:
            for family in combined:
                for table in family:
                    table.close()
            for table in rounded.values():
                table.close()
        for a in (kinetic,electron,photon,transport,electrons,emission):
            _readonly(a)
        return SolverResult(kinetic,electron,photon,
            {name:transport[:,i,:] for i,name in enumerate(TRANSPORT)},
            {name:electrons[:,i,:] for i,name in enumerate(ELECTRONS)},
            {name:emission[:,i,:] for i,name in enumerate(EMISSION)},
            {"solver_abi": 2, "legacy_table_precision": legacy_table_precision,
             "catalogue": rows.tolist(), "catalogue_columns": ["z","Mstar_Msun","Re_kpc","SFR_Msun_yr"],
             "cosmic_ray_bins": ne, "photon_bins": np_, "solver_cells": grid.cells,
             "table_grid": [preparation.grid.nx, preparation.grid.ny],
             "threads": threads, "preparation_fingerprint": preparation._library_hash,
             "solver_sha256": hashlib.sha256(Path(lib._name).read_bytes()).hexdigest(),
             "emission_units": "GeV^-1 s^-1 except dimensionless tau_free_free",
             "particle_units": "GeV^-1", "frame": "source; disc free-free applied"})
