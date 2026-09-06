"""Week-2 preparation: orchestration and file I/O in Python, physics in C."""
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
from contextlib import ExitStack
import numpy as np

from ._bindings import DoublePointer as DP, check
from .inputs import Galaxy, Grid
from .model import Context

P = ct.c_void_p
SP = ct.POINTER(ct.c_size_t)
PP = ct.POINTER(P)
FIELDS = ("3000", "4000", "7500", "UV", "CMB", "FIR")
KINDS = ("emission", "gamma", "bs", "sy")
PROPERTY_NAMES = ("height_pc", "density_cm3", "magnetic_field_gauss",
                  "gas_dispersion_kms", "area_pc2", "gas_surface_msun_pc2",
                  "sfr_surface_msun_yr_pc2", "stellar_surface_msun_pc2",
                  "dust_temperature_k", "halo_magnetic_field_gauss")
RADIATION_NAMES = ("CMB", "FIR", "3000", "4000", "7500", "UV", "total")

def bind(lib):
    signatures = {
        "cg_galaxies": ([P, ct.c_size_t, DP, DP], ct.c_int),
        "cg_preparation_bounds": ([DP, ct.c_size_t, DP, DP], ct.c_int),
        "cg_radiation": ([DP, ct.c_size_t, DP, DP, DP, DP], ct.c_int),
        "cg_table_create": ([ct.c_int, ct.c_int, ct.c_double, ct.c_size_t,
                             ct.c_size_t, DP, PP], ct.c_int),
        "cg_table_import": ([ct.c_size_t, ct.c_size_t, DP, DP, DP, PP], ct.c_int),
        "cg_table_borrow": ([ct.c_size_t, ct.c_size_t, DP, DP, DP, PP], ct.c_int),
        "cg_table_destroy": ([P], None),
        "cg_table_shape": ([P, SP, SP], ct.c_int),
        "cg_table_copy": ([P, DP, DP, DP], ct.c_int),
        "cg_table_eval": ([P, ct.c_size_t, DP, DP, DP], ct.c_int),
        "cg_temperature_grid": ([ct.c_int, ct.c_double, ct.c_double, ct.c_double,
                                  ct.c_size_t, DP, SP], ct.c_int),
        "cg_combine_ic": ([DP, PP, ct.c_double, ct.c_double, PP], ct.c_int),
    }
    for name, (args, result) in signatures.items():
        f = getattr(lib, name)
        f.argtypes, f.restype = args, result

def doubles(values):
    values = tuple(values)
    return (ct.c_double * len(values))(*values)

@dataclass(frozen=True)
class TableData:
    x: tuple
    y: tuple
    values: tuple  # flattened [iy*nx+ix], never transposed

class Table:
    """Python-owned NumPy storage, borrowed by a small C evaluation descriptor.

    Construction consumes a temporary C-owned handle; only the unchanged
    legacy generation scratch is C-owned. Persistent arrays belong to NumPy.
    """
    def __init__(self, lib, handle, metadata):
        self._lib, self._handle = lib, P()
        self.metadata = dict(metadata)
        self._lock = threading.RLock()
        try:
            nx, ny = ct.c_size_t(), ct.c_size_t()
            check(lib, lib.cg_table_shape(handle, ct.byref(nx), ct.byref(ny)))
            self._x = np.empty(nx.value, dtype=np.float64)
            self._y = np.empty(ny.value, dtype=np.float64)
            self._values = np.empty((ny.value, nx.value), dtype=np.float64)
            check(lib, lib.cg_table_copy(handle, self._pointer(self._x),
                  self._pointer(self._y), self._pointer(self._values)))
            self._borrow()
        finally:
            lib.cg_table_destroy(handle)

    @staticmethod
    def _pointer(values):
        return values.ctypes.data_as(DP)

    def _borrow(self):
        # Strong references on self outlive the C descriptor. Read-only
        # arrays prevent accidental edits while a C call has released the GIL.
        for values in (self._x, self._y, self._values):
            values.flags.writeable = False
        check(self._lib, self._lib.cg_table_borrow(
            self._x.size, self._y.size, self._pointer(self._x),
            self._pointer(self._y), self._pointer(self._values),
            ct.byref(self._handle)))

    @classmethod
    def from_data(cls, lib, data, metadata):
        """Load Python cache data directly into NumPy, without C-owned copies."""
        table = cls.__new__(cls)
        table._lib, table._handle = lib, P()
        table._lock = threading.RLock()
        table.metadata = dict(metadata)
        table._x = np.array(data.x, dtype=np.float64, order="C", copy=True)
        table._y = np.array(data.y, dtype=np.float64, order="C", copy=True)
        if table._x.ndim != 1 or table._y.ndim != 1:
            raise ValueError("Table axes must be one-dimensional")
        values = np.asarray(data.values, dtype=np.float64)
        if values.size != table._x.size * table._y.size:
            raise ValueError("Table values must contain nx * ny entries")
        table._values = np.empty((table._y.size, table._x.size), dtype=np.float64)
        table._values[:] = values.reshape(table._values.shape)
        table._borrow()
        return table

    def _view(self, values):
        with self._lock:
            self._ensure_open()
            return values.view()

    @property
    def x(self):
        """Read-only NumPy axis in the units recorded in metadata."""
        return self._view(self._x)

    @property
    def y(self):
        return self._view(self._y)

    @property
    def values(self):
        """Read-only C-contiguous NumPy array shaped (ny, nx)."""
        return self._view(self._values)

    def _ensure_open(self):
        if not self._handle.value:
            raise RuntimeError("Table is closed")

    def snapshot(self):
        with self._lock:
            self._ensure_open()
            return TableData(tuple(self._x), tuple(self._y), tuple(self._values.ravel()))

    def evaluate(self, x, y=None):
        with self._lock:
            self._ensure_open()
            xs = doubles(x)
            ys = None if y is None else doubles(y)
            if ys is not None and len(xs) != len(ys):
                raise ValueError("x and y lengths must match")
            result = doubles([0.] * len(xs))
            check(self._lib, self._lib.cg_table_eval(self._handle, len(xs), xs, ys, result))
            return list(result)

    def close(self):
        with self._lock:
            if self._handle.value:
                self._lib.cg_table_destroy(self._handle)
                self._handle = P()

    def __enter__(self):
        self._ensure_open()
        return self

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        if getattr(self, "_handle", None):
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
    files. All temperature interpolation/physics arithmetic stays in C.
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
        self._library_hash = hashlib.sha256(Path(self._lib._name).read_bytes()).hexdigest()
        rows = doubles(v for g in self.galaxies for v in g.as_row())
        output = doubles([0.] * (10 * len(self.galaxies)))
        check(self._lib, self._lib.cg_galaxies(self._context._handle, len(self.galaxies), rows, output))
        self.properties = tuple(dict(zip(PROPERTY_NAMES, output[i:i+10]))
                                for i in range(0, len(output), 10))
        config, temperatures = doubles([0.] * 8), doubles([0.] * 3)
        check(self._lib, self._lib.cg_preparation_bounds(rows, len(self.galaxies), config, temperatures))
        self.config = tuple(config)
        self.temperature_bounds = tuple(temperatures)

    def _ensure_open(self):
        if not self._context._handle.value:
            raise RuntimeError("Preparation is closed")

    def temperature_grid(self, field):
        with self._lock:
            self._ensure_open()
            if field == "CMB":
                args = (4, 0., self.temperature_bounds[0], .5)
            elif field == "FIR":
                args = (5, self.temperature_bounds[1], self.temperature_bounds[2], 5.)
            else:
                raise ValueError("Temperature grid is only for CMB or FIR")
            count = ct.c_size_t()
            check(self._lib, self._lib.cg_temperature_grid(*args, 0, None, ct.byref(count)))
            values = doubles([0.] * count.value)
            check(self._lib, self._lib.cg_temperature_grid(*args, count.value, values, ct.byref(count)))
            return tuple(values)

    def radiation(self, index, energies):
        with self._lock:
            self._ensure_open()
            g = self.galaxies[index]
            es = doubles(energies)
            fields, urad = doubles([0.] * (len(es) * 7)), doubles([0.] * 8)
            check(self._lib, self._lib.cg_radiation(doubles(g.as_row()), len(es), es,
                  doubles(self.config[4:6]), fields, urad))
            return {"fields_cm3_gev": {name: tuple(fields[j::7])
                    for j, name in enumerate(RADIATION_NAMES)}, "urad_ub_ev_cm3": tuple(urad)}

    def table(self, kind, field="3000", temperature=0.):
        with self._lock:
            self._ensure_open()
            if kind not in KINDS or field not in FIELDS:
                raise ValueError("Unknown table kind/field")
            if kind in ("bs", "sy"):
                field, temperature = "3000", 0.
            if field not in ("CMB", "FIR"):
                temperature = 0.
            meta = {"schema": 1, "library_sha256": self._library_hash,
                    "kind": kind, "field": field, "temperature_k": float(temperature),
                    "nx": self.grid.nx, "ny": self.grid.ny, "config": list(self.config),
                    "units": "dimensionless" if kind == "sy" else "mb/GeV" if kind == "bs" else "1/(s GeV)",
                    "x_axis": "x" if kind == "sy" else "DeltaE_GeV" if kind == "gamma" else "photon_GeV",
                    "y_axis": "unused" if kind == "sy" else "electron_total_GeV"}
            encoded = json.dumps(meta, sort_keys=True, allow_nan=False).encode()
            key = hashlib.sha256(encoded).hexdigest()
            if key in self._tables:
                return self._tables[key]
            handle = P()
            path = None if self.cache is None else self.cache / (key + ".cgt")
            if path is not None and path.exists():
                data = _read_cache(path, meta)
                table = Table.from_data(self._lib, data, meta)
                self.cache_hits += 1
            else:
                check(self._lib, self._lib.cg_table_create(KINDS.index(kind), FIELDS.index(field),
                      temperature, self.grid.nx, self.grid.ny, doubles(self.config), ct.byref(handle)))
                table = Table(self._lib, handle, meta)
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
        """Independent Python-owned arrays with a borrowed C evaluation handle."""
        with self._lock:
            self._ensure_open()
            if kind not in ("emission", "gamma"):
                raise ValueError("Expected emission or gamma")
            galaxy = self.galaxies[index]
            planes = [self.table(kind, field) for field in FIELDS[:4]]
            fractions = []
            # Native C computes both CMB temperature and dust temperature.
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
            handle = P()
            with ExitStack() as locks:
                for table in planes:
                    locks.enter_context(table._lock)
                    table._ensure_open()
                check(self._lib, self._lib.cg_combine_ic(doubles(galaxy.as_row()),
                      (P * 8)(*(t._handle for t in planes)), *fractions, ct.byref(handle)))
            return Table(self._lib, handle, {"kind": kind, "combination": "legacy-reversed-temperature-weights",
                         "units": "1/(s GeV)", "x_axis": "DeltaE_GeV" if kind == "gamma" else "photon_GeV"})

    def _single_cmb_temperature(self, galaxy):
        config, temps = doubles([0.] * 8), doubles([0.] * 3)
        check(self._lib, self._lib.cg_preparation_bounds(doubles(galaxy.as_row()), 1, config, temps))
        return temps[0]

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
