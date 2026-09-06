"""Week-1 context only; not yet the end-to-end galaxy model."""
import ctypes as ct
import threading
from ._bindings import load_library, check

class Context:
    """C-owned configuration; use as a context manager.

    ionisation_loss takes an iterable of TOTAL electron energies in GeV
    and a hydrogen number density in cm^-3. Returns a list in GeV/s.
    Inputs are copied to contiguous doubles; C retains no Python pointers.
    No physics constants or equations are reimplemented in Python.
    """
    def __init__(self, threads=1, library=None):
        if isinstance(threads, bool) or not isinstance(threads, int):
            raise TypeError("threads must be an integer")
        if not 1 <= threads <= 2147483647:
            raise ValueError("threads must fit a positive C int")
        self._lock = threading.RLock()
        self._handle = ct.c_void_p()
        self._lib = load_library(library)
        check(self._lib, self._lib.cg_context_create(threads, ct.byref(self._handle)))

    @property
    def version(self):
        return self._lib.cg_version().decode()

    @property
    def openmp_enabled(self):
        return bool(self._lib.cg_openmp_enabled())

    def ionisation_loss(self, energies_gev, density_cm3):
        with self._lock:
            if not self._handle.value:
                raise RuntimeError("Context is closed")
            values = tuple(float(value) for value in energies_gev)
            buffer_type = ct.c_double * len(values)
            energies = buffer_type(*values)
            output = buffer_type()
            check(self._lib, self._lib.cg_ionisation(
                self._handle, len(values), energies, float(density_cm3), output))
            return list(output)

    def close(self):
        with self._lock:
            if self._handle.value:
                self._lib.cg_context_destroy(self._handle)
                self._handle = ct.c_void_p()

    def __enter__(self):
        if not self._handle.value:
            raise RuntimeError("Context is closed")
        return self

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        if getattr(self, "_handle", None) and hasattr(self, "_lib"):
            self.close()

