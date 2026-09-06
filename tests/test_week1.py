import ctypes as ct
import math
from pathlib import Path
import subprocess
import unittest
from concurrent.futures import ThreadPoolExecutor
from congruents import Context
from congruents._bindings import load_library

ROOT = Path(__file__).resolve().parents[1]

class Week1Tests(unittest.TestCase):
    def test_direct_c_and_threads(self):
        reference = [float(x) for x in subprocess.check_output(
            [str(ROOT / "build/direct_ionisation")], text=True).split()]
        for threads in (1, 4):
            with Context(threads) as context:
                self.assertEqual(context.version, "0.1.0")
                self.assertTrue(context.openmp_enabled)
                actual = []
                for density in (0., 1.e-3, 1., 1.e3):
                    actual.extend(context.ionisation_loss([.001, .01, 1., 100., 1.e5], density))
                for got, expected in zip(actual, reference):
                    self.assertTrue(math.isclose(got, expected, rel_tol=1e-14, abs_tol=0))
                self.assertEqual(len(actual), len(reference))

    def test_lifecycle(self):
        for _ in range(1000):
            with Context() as context:
                self.assertEqual(context.ionisation_loss([], 1.), [])
            context.close()
            with self.assertRaises(RuntimeError):
                context.ionisation_loss([1.], 1.)

    def test_validation_survives(self):
        with Context() as context:
            for energies, density in [([0.], 1.), ([.0001], 1.),
                ([float("nan")], 1.), ([float("inf")], 1.),
                ([1.], -1.), ([1.], float("nan"))]:
                with self.assertRaises(ValueError):
                    context.ionisation_loss(energies, density)
            self.assertLess(context.ionisation_loss([1.], 1.)[0], 0)
        for threads in (0, -1, 2**40):
            with self.assertRaises(ValueError):
                Context(threads)
        with self.assertRaises(TypeError):
            Context(1.5)

    def test_missing_library(self):
        with self.assertRaises(FileNotFoundError):
            Context(library=ROOT / "does-not-exist.so")

    def test_numerical_error_recovery(self):
        with Context() as context:
            with self.assertRaises(RuntimeError):
                context.ionisation_loss([1.e308], 1.)
            self.assertLess(context.ionisation_loss([1.], 1.)[0], 0)

    def test_concurrent_contexts(self):
        def evaluate(threads):
            with Context(threads) as context:
                return context.ionisation_loss([.01, 1., 100.] * 100, 10.)
        expected = evaluate(1)
        with ThreadPoolExecutor(max_workers=4) as pool:
            for actual in pool.map(evaluate, [1, 2, 3, 4] * 4):
                self.assertEqual(actual, expected)

    def test_raw_boundary(self):
        lib = load_library()
        handle = ct.c_void_p()
        self.assertEqual(lib.cg_context_create(0, ct.byref(handle)), 1)
        self.assertFalse(handle.value)
        self.assertEqual(lib.cg_context_create(1, None), 1)
        self.assertEqual(lib.cg_context_create(1, ct.byref(handle)), 0)
        try:
            out = (ct.c_double * 1)(123.)
            bad = (ct.c_double * 1)(float("nan"))
            self.assertEqual(lib.cg_ionisation(handle, 1, bad, 1., out), 1)
            self.assertEqual(out[0], 123.)
            self.assertEqual(lib.cg_ionisation(handle, 1, None, 1., out), 1)
            self.assertEqual(lib.cg_ionisation(None, 0, None, 1., None), 1)
        finally:
            lib.cg_context_destroy(handle)
            lib.cg_context_destroy(None)

if __name__ == "__main__":
    unittest.main()
