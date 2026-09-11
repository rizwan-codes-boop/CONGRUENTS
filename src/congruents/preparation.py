"""Python serial preparation with a native, galaxy-parallel property loop."""
from array import array
import bisect
import ctypes as ct
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import tempfile
import threading
from . import serial, constants
import scipy
import numpy as np

from ._bindings import DoublePointer as DP, check
from .inputs import Galaxy, Grid
from .model import Context

P = ct.c_void_p
FIELDS = ("3000", "4000", "7500", "UV", "CMB", "FIR")
KINDS = ("emission", "gamma", "bs", "sy")
PROPERTY_NAMES = ("height_pc", "density_cm3", "magnetic_field_gauss",
                  "gas_dispersion_kms", "area_pc2", "gas_surface_msun_pc2",
                  "sfr_surface_msun_yr_pc2", "stellar_surface_msun_pc2",
                  "dust_temperature_k", "halo_magnetic_field_gauss")
RADIATION_NAMES = ("CMB", "FIR", "3000", "4000", "7500", "UV", "total")

def bind(lib):
    lib.cg_galaxies.argtypes = [P, ct.c_size_t, DP, DP]
    lib.cg_galaxies.restype = ct.c_int

def doubles(values):
    values = tuple(values)
    return (ct.c_double * len(values))(*values)

@dataclass(frozen=True)
class TableData:
    x: tuple
    y: tuple
    values: tuple  # flattened [iy*nx+ix], never transposed

class Table:
    """Read-only Python/NumPy tables; serial interpolation stays in Python."""
    def __init__(self, data, metadata):
        self.metadata = dict(metadata)
        self._lock = threading.RLock()
        self._closed = False
        self._x = np.array(data.x, dtype=np.float64, copy=True)
        self._y = np.array(data.y, dtype=np.float64, copy=True)
        if self._x.ndim != 1 or self._y.ndim != 1:
            raise ValueError("Table axes must be one-dimensional")
        nx, ny = self._x.size, self._y.size
        if not 2 <= nx <= 4096 or not 1 <= ny <= 4096 or nx*ny > 4194304:
            raise ValueError("Invalid table shape")
        for a in (self._x,self._y):
            if not np.isfinite(a).all() or (a<=0).any() or (np.diff(a)<=0).any():
                raise ValueError("Table axes must be finite, positive and increasing")
        values = np.asarray(data.values, dtype=np.float64)
        if values.size != nx*ny or not np.isfinite(values).all() or (values<0).any():
            raise ValueError("Invalid table values")
        self._values = np.empty((ny,nx),dtype=np.float64)
        self._values[:] = values.reshape(ny,nx)
        for a in (self._x,self._y,self._values):
            a.flags.writeable = False

    @classmethod
    def from_data(cls, lib, data, metadata):
        # Retained Python compatibility; no C handle is created.
        return cls(data,metadata)

    def _ensure_open(self):
        if self._closed:
            raise RuntimeError("Table is closed")

    def _view(self, a):
        with self._lock:
            self._ensure_open()
            return a.view()

    @property
    def x(self):
        return self._view(self._x)

    @property
    def y(self):
        return self._view(self._y)

    @property
    def values(self):
        return self._view(self._values)

    def snapshot(self):
        with self._lock:
            self._ensure_open()
            return TableData(tuple(self._x),tuple(self._y),tuple(self._values.ravel()))

    def evaluate(self, x, y=None):
        with self._lock:
            self._ensure_open()
            x=np.asarray(tuple(x),dtype=np.float64)
            if x.ndim!=1 or not np.isfinite(x).all():
                raise ValueError("Coordinates must be finite one-dimensional arrays")
            if len(self._y)==1:
                return np.interp(x,self._x,self._values[0],left=0,right=0).tolist()
            if y is None:
                raise ValueError("A two-dimensional table needs electron coordinates")
            y=np.asarray(tuple(y),dtype=np.float64)
            if y.shape!=x.shape or not np.isfinite(y).all():
                raise ValueError("Coordinate lengths must match and be finite")
            valid=(x>=self._x[0])&(x<=self._x[-1])&(y>=self._y[0])&(y<=self._y[-1])
            out=np.zeros_like(x)
            xv,yv=x[valid],y[valid]
            ix=np.clip(np.searchsorted(self._x,xv,side="right")-1,0,len(self._x)-2)
            iy=np.clip(np.searchsorted(self._y,yv,side="right")-1,0,len(self._y)-2)
            u=(xv-self._x[ix])/(self._x[ix+1]-self._x[ix])
            v=(yv-self._y[iy])/(self._y[iy+1]-self._y[iy])
            z=self._values
            out[valid]=(1-v)*((1-u)*z[iy,ix]+u*z[iy,ix+1])+v*((1-u)*z[iy+1,ix]+u*z[iy+1,ix+1])
            return out.tolist()

    def close(self):
        with self._lock:
            self._closed=True

    def __enter__(self):
        self._ensure_open()
        return self

    def __exit__(self,*args):
        self.close()

def _save_cache(path, metadata, data):
    payload = array("d", data.x + data.y + data.values)
    if sys.byteorder != "little":
        payload.byteswap()
    raw = payload.tobytes()
    header = json.dumps({"metadata": metadata, "nx": len(data.x), "ny": len(data.y),
                         "sha256": hashlib.sha256(raw).hexdigest()}, sort_keys=True).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic replacement; an interrupted writer cannot leave a partial cache.
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            stream.write(b"CGT1" + struct.pack("<I", len(header)) + header + raw)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    try:
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)

def _read_cache(path, expected):
    with path.open("rb") as stream:
        prefix = stream.read(8)
        if len(prefix) != 8 or prefix[:4] != b"CGT1":
            raise ValueError("Invalid cache header")
        length = struct.unpack("<I", prefix[4:])[0]
        if length > 65536:
            raise ValueError("Oversized cache header")
        header = json.loads(stream.read(length))
        if header["metadata"] != expected:
            raise ValueError("Cache metadata mismatch")
        nx, ny = header["nx"], header["ny"]
        expected_ny = 1 if expected["kind"] == "sy" else expected["ny"]
        if nx != expected["nx"] or ny != expected_ny:
            raise ValueError("Cache shape mismatch")
        size = (nx + ny + nx * ny) * 8
        raw = stream.read(size + 1)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != header["sha256"]:
            raise ValueError("Cache payload length/checksum mismatch")
    values = array("d")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    return TableData(tuple(values[:nx]), tuple(values[nx:nx+ny]), tuple(values[nx+ny:]))

class Preparation:
    """Own reusable per-field tables for a validated catalogue.

    Defaults reproduce the production grid sizes; use Grid(8,8) ONLY for smoke
    tests. Cache files are new versioned binary artifacts, never legacy text
    files. Serial generation and evaluation are Python; galaxy batches use C.
    """
    def __init__(self, galaxies, grid=None, threads=1, cache=None, library=None):
        self.galaxies = tuple(galaxies)
        if not self.galaxies or len(self.galaxies) > 100000:
            raise ValueError("Require 1..100000 galaxies")
        if not all(isinstance(g, Galaxy) for g in self.galaxies):
            raise TypeError("Expected Galaxy objects")
        self.grid = grid if grid is not None else Grid()
        if not isinstance(self.grid, Grid):
            raise TypeError("grid must be Grid")
        self._lock = threading.RLock()
        self._context = Context(threads, library)
        self._lib = self._context._lib
        bind(self._lib)
        self._tables = {}
        self.cache = None if cache is None else Path(cache)
        self.cache_hits = self.cache_misses = 0
        source = b"".join(Path(__file__).with_name(name).read_bytes() for name in
                          ("serial.py","constants.py","quadrature.py","preparation.py"))
        self._library_hash = hashlib.sha256(Path(self._lib._name).read_bytes()+source+
                          (np.__version__+scipy.__version__).encode()).hexdigest()
        rows = doubles(v for g in self.galaxies for v in g.as_row())
        output = doubles([0.] * (10 * len(self.galaxies)))
        check(self._lib, self._lib.cg_galaxies(self._context._handle, len(self.galaxies), rows, output))
        self.properties = tuple(dict(zip(PROPERTY_NAMES, output[i:i+10]))
                                for i in range(0, len(output), 10))
        self.config, self.temperature_bounds = serial.bounds(self.galaxies,self.properties)

    def _ensure_open(self):
        if not self._context._handle.value:
            raise RuntimeError("Preparation is closed")

    def temperature_grid(self, field):
        with self._lock:
            self._ensure_open()
            return serial.temperature_grid(field,self.temperature_bounds)

    def radiation(self,index,energies):
        with self._lock:
            self._ensure_open()
            fields,urad=serial.radiation(self.galaxies[index],self.properties[index],energies,self.config)
            return {"fields_cm3_gev":dict(zip(RADIATION_NAMES,fields)),"urad_ub_ev_cm3":urad}

    def table(self, kind, field="3000", temperature=0.):
        with self._lock:
            self._ensure_open()
            if kind not in KINDS or field not in FIELDS:
                raise ValueError("Unknown table kind/field")
            if kind in ("bs", "sy"):
                field, temperature = "3000", 0.
            if field not in ("CMB", "FIR"):
                temperature = 0.
            meta = {"schema": 2, "implementation_sha256": self._library_hash,
                    "kind": kind, "field": field, "temperature_k": float(temperature),
                    "nx": self.grid.nx, "ny": self.grid.ny, "config": list(self.config),
                    "units": "dimensionless" if kind == "sy" else "mb/GeV" if kind == "bs" else "1/(s GeV)",
                    "x_axis": "x" if kind == "sy" else "DeltaE_GeV" if kind == "gamma" else "photon_GeV",
                    "y_axis": "unused" if kind == "sy" else "electron_total_GeV"}
            encoded = json.dumps(meta, sort_keys=True, allow_nan=False).encode()
            key = hashlib.sha256(encoded).hexdigest()
            if key in self._tables:
                return self._tables[key]
            path = None if self.cache is None else self.cache / (key + ".cgt")
            if path is not None and path.exists():
                data = _read_cache(path, meta)
                table = Table.from_data(self._lib, data, meta)
                self.cache_hits += 1
            else:
                x,y,z=serial.generate(kind,field,temperature,self.grid.nx,self.grid.ny,self.config)
                table=Table(TableData(x,y,z),meta)
                self.cache_misses += 1
            try:
                if path is not None and not path.exists():
                    _save_cache(path, meta, table.snapshot())
            except BaseException:
                table.close()
                raise
            self._tables[key] = table
            return table

    def prepare_tables(self):
        """Prepare all 14 logical table families; CMB/FIR contain multiple planes."""
        with self._lock:
            for kind in ("emission", "gamma"):
                for field in FIELDS[:4]:
                    self.table(kind, field)
                for field in ("CMB", "FIR"):
                    for temperature in self.temperature_grid(field):
                        self.table(kind, field, temperature)
            self.table("bs")
            self.table("sy")
            return len(self._tables)

    def combined_ic(self, index, kind="emission"):
        """Independent Python-combined table; legacy reversed weights are retained."""
        with self._lock:
            self._ensure_open()
            if kind not in ("emission", "gamma"):
                raise ValueError("Expected emission or gamma")
            galaxy = self.galaxies[index]
            planes = [self.table(kind, field) for field in FIELDS[:4]]
            fractions = []
            # C computed dust temperature in its galaxy loop; CMB setup is serial Python.
            target_cmb = self._single_cmb_temperature(galaxy)
            for field, target in (("CMB", target_cmb),
                                  ("FIR", self.properties[index]["dust_temperature_k"])):
                grid = self.temperature_grid(field)
                j = bisect.bisect_right(grid, target) - 1
                if j < 0 or j >= len(grid)-1:
                    raise ValueError("Temperature has no safe legacy interpolation bracket")
                # Match legacy position-index interpolation before taking fraction.
                position = j + (target-grid[j])/(grid[j+1]-grid[j])
                fractions.append(position-j)
                planes.extend([self.table(kind, field, grid[j]), self.table(kind, field, grid[j+1])])
            for table in planes:
                table._ensure_open()
            d=serial.dilution(galaxy,self.properties[index]["dust_temperature_k"])
            fc,ff=fractions
            z=[table.values for table in planes]
            result=(z[4]*fc+z[5]*(1-fc))+d[4]*(z[6]*ff+z[7]*(1-ff))
            result=result+d[0]*z[0]+d[1]*z[1]+d[2]*z[2]+d[3]*z[3]
            return Table(TableData(planes[0].x,planes[0].y,result),
                         {"kind":kind,"combination":"legacy-reversed-temperature-weights",
                          "units":"1/(s GeV)","x_axis":"DeltaE_GeV" if kind=="gamma" else "photon_GeV"})

    def _single_cmb_temperature(self,galaxy):
        return constants.TCMB*(1+galaxy.redshift)

    def close(self):
        with self._lock:
            for table in self._tables.values():
                table.close()
            self._tables.clear()
            self._context.close()

    def __enter__(self):
        self._ensure_open()
        return self

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        if hasattr(self, "_tables"):
            self.close()
