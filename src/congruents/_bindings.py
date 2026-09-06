"""Explicit ctypes ABI, following the SLUG bayesphot loading pattern."""
import ctypes as ct
import os
from pathlib import Path
import sys

DoublePointer = ct.POINTER(ct.c_double)

def load_library(path=None):
    if path is None:
        path = os.environ.get("CONGRUENTS_LIBRARY")
    if path is None:
        suffix = ".dylib" if sys.platform == "darwin" else ".so"
        path = Path(__file__).resolve().parents[2] / "build" / ("libcongruents" + suffix)
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"C library not found: {path}. Run make shared or set CONGRUENTS_LIBRARY."
        )
    lib = ct.CDLL(str(path))
    signatures = {
        "cg_version": (ct.c_char_p, []),
        "cg_abi_version": (ct.c_uint, []),
        "cg_status_message": (ct.c_char_p, [ct.c_int]),
        "cg_openmp_enabled": (ct.c_int, []),
        "cg_context_create": (ct.c_int, [ct.c_int, ct.POINTER(ct.c_void_p)]),
        "cg_context_destroy": (None, [ct.c_void_p]),
        "cg_ionisation": (ct.c_int, [ct.c_void_p, ct.c_size_t,
                                    DoublePointer, ct.c_double, DoublePointer]),
    }
    for name, (result, args) in signatures.items():
        function = getattr(lib, name)
        function.restype = result
        function.argtypes = args
    if lib.cg_abi_version() != 1:
        raise RuntimeError("Unsupported CONGRUENTS C ABI")
    return lib

def check(lib, status):
    if status:
        message = lib.cg_status_message(status).decode()
        if status == 1:
            raise ValueError(message)
        if status == 2:
            raise MemoryError(message)
        raise RuntimeError(message)

