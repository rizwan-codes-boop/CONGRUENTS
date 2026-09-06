import ctypes as ct
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
import gc
import numpy as np

from congruents import Galaxy, Grid, Preparation, load_catalogue, write_catalogue
from congruents.preparation import KINDS, FIELDS, PROPERTY_NAMES, Table, TableData, doubles
from congruents._bindings import check

ROOT = Path(__file__).resolve().parents[1]
GALAXIES = load_catalogue(ROOT / "input/cat_nt.txt")

class Week2Tests(unittest.TestCase):
    def assertClose(self, actual, expected, tolerance=1e-12):
        self.assertEqual(len(actual), len(expected))
        for i, (a, b) in enumerate(zip(actual, expected)):
            self.assertTrue(math.isclose(a, b, rel_tol=tolerance, abs_tol=1e-280),
                            (i, a, b))

    def reference(self, p, kind, field, temperature=0., galaxy=None):
        args = [str(ROOT/"build/direct_preparation"), str(KINDS.index(kind)),
                str(6 if galaxy is not None else FIELDS.index(field)), str(temperature),
                str(p.grid.nx), str(p.grid.ny), *(repr(v) for v in p.config)]
        if galaxy is not None:
            args.extend(repr(v) for v in galaxy.as_row()+p.temperature_bounds)
        words = subprocess.check_output(args, text=True).split()
        nx, ny = map(int, words[:2])
        vals = tuple(map(float, words[2:]))
        self.assertEqual(len(vals), nx+ny+nx*ny)
        return TableData(vals[:nx], vals[nx:nx+ny], vals[nx+ny:])

    def test_catalogue_roundtrip_and_validation(self):
        self.assertEqual(len(GALAXIES), 11)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"cat.txt"
            write_catalogue(path, GALAXIES)
            self.assertEqual(load_catalogue(path), GALAXIES)
            path.write_text("n_gal\n2\nz Mstar__Msol Re__kpc SFR__Msolyrm1\n0 1 1 1\n")
            with self.assertRaises(ValueError):
                load_catalogue(path)
        for row in ((0, 1, 0, 1), (0, 1, 1, 0), (-1, 1, 1, 1), (0, float("nan"), 1, 1)):
            with self.assertRaises(ValueError):
                Galaxy(*row)
        for n in (0, 1, 4097):
            with self.assertRaises(ValueError):
                Grid(n, 8)
        with self.assertRaises(TypeError):
            Grid(2.5, 8)

    @unittest.skipUnless(importlib.util.find_spec("astropy"), "optional Astropy extra")
    def test_astropy_units(self):
        from astropy import units as u
        g = Galaxy.from_quantities(.01, 1e10*u.Msun, 2000*u.pc, 2*u.Msun/u.yr)
        self.assertEqual(g.radius_kpc, 2.)
        with self.assertRaises(u.UnitConversionError):
            Galaxy.from_quantities(.01, 1e10*u.Msun, 2*u.s, 2*u.Msun/u.yr)

    def test_properties_and_radiation_baseline(self):
        def fixture(name):
            return [list(map(float, row.split())) for row in
                    (ROOT/"tests/fixtures"/name).read_text().splitlines()[1:]]
        expected_props, expected_rad = fixture("gal_data.txt"), fixture("Urad_Ub.txt")
        for threads in (1, 4):
            with Preparation(GALAXIES, Grid(8, 8), threads=threads) as p:
                for i, prop in enumerate(p.properties):
                    self.assertClose([prop[k] for k in PROPERTY_NAMES[:9]], expected_props[i], 6e-7)
                    result = p.radiation(i, [1e-12, 1e-9])
                    self.assertClose(result["urad_ub_ev_cm3"], expected_rad[i], 6e-7)
                    self.assertTrue(all(v >= 0 and math.isfinite(v)
                        for field in result["fields_cm3_gev"].values() for v in field))

    def test_all_table_families_direct_c(self):
        with Preparation(GALAXIES, Grid(8, 8)) as p:
            for kind in ("emission", "gamma"):
                for field in FIELDS:
                    temps = p.temperature_grid(field) if field in ("CMB", "FIR") else (0.,)
                    for temp in temps:
                        with self.subTest(kind=kind, field=field, temperature=temp):
                            actual = p.table(kind, field, temp).snapshot()
                            expected = self.reference(p, kind, field, temp)
                            self.assertEqual(actual, expected)
                            self.assertEqual(len(actual.values), 64)
                            self.assertTrue(all(x >= 0 and math.isfinite(x) for x in actual.values))
            for kind in ("bs", "sy"):
                self.assertEqual(p.table(kind).snapshot(), self.reference(p, kind, "3000"))
            self.assertEqual(p.prepare_tables(), 26)
            emission = p.table("emission").snapshot()
            gamma = p.table("gamma").snapshot()
            self.assertTrue(math.isclose(emission.x[0], p.config[0], rel_tol=1e-14))
            self.assertTrue(math.isclose(gamma.x[0], p.config[4], rel_tol=1e-14))
            self.assertNotEqual(emission.values, gamma.values)

    def test_combined_ic_direct_legacy(self):
        with Preparation(GALAXIES, Grid(8, 8)) as p:
            for kind in ("emission", "gamma"):
                for i in (0, 1, 7, 10):
                    with p.combined_ic(i, kind) as table:
                        got = table.snapshot()
                        expected = self.reference(p, kind, "3000", galaxy=GALAXIES[i])
                        self.assertEqual(got.x, expected.x)
                        self.assertEqual(got.y, expected.y)
                        self.assertClose(got.values, expected.values, 2e-13)
                        # Test interpolation nodes and outside-domain zero.
                        self.assertClose(table.evaluate([got.x[3]], [got.y[3]]),
                                         [got.values[3*len(got.x)+3]])
                        self.assertEqual(table.evaluate([0.], [1.]), [0.])

    def test_production_resolution_spot(self):
        fixture = json.loads((ROOT/"tests/fixtures/production_ic_spot.json").read_text())
        with Preparation(GALAXIES, Grid(2, 2)) as p:
            config = list(p.config)
            config[:2], config[2:4] = fixture["x"], fixture["y"]
            handle = ct.c_void_p()
            check(p._lib, p._lib.cg_table_create(0, 0, 0., 2, 2,
                  doubles(config), ct.byref(handle)))
            with Table(p._lib, handle, {}) as table:
                # Saved axes/values have only 7 significant digits.
                self.assertClose(table.snapshot().values, fixture["values"], 1e-5)

    def test_cache_cold_warm_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                p.prepare_tables()
                self.assertEqual((p.cache_hits, p.cache_misses), (0, 26))
                reference = {k: t.snapshot() for k, t in p._tables.items()}
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                p.prepare_tables()
                self.assertEqual((p.cache_hits, p.cache_misses), (26, 0))
                self.assertEqual({k: t.snapshot() for k, t in p._tables.items()}, reference)
            path = next(Path(tmp).glob("*.cgt"))
            raw = bytearray(path.read_bytes())
            raw[-1] ^= 1
            path.write_bytes(raw)
            with Preparation(GALAXIES, Grid(8, 8), cache=tmp) as p:
                with self.assertRaises(ValueError):
                    p.prepare_tables()
            # Grid changes produce distinct keys, not reuse of wrong-shaped tables.
            with Preparation(GALAXIES, Grid(9, 8), cache=tmp) as p:
                p.table("bs")
                self.assertEqual(p.cache_misses, 1)

    def test_lifetime_and_invalid_inputs(self):
        p = Preparation(GALAXIES, Grid(8, 8))
        with self.assertRaises(ValueError):
            p.radiation(0, [float("nan")])
        with self.assertRaises(ValueError):
            p.table("emission", "CMB", -1)
        t = p.table("bs")
        p.close()
        p.close()
        with self.assertRaises(RuntimeError):
            t.snapshot()
        with self.assertRaises(RuntimeError):
            p.radiation(0, [1e-9])
        # Legacy single-plane temperatures are rejected instead of underflowing.
        with Preparation([GALAXIES[0]], Grid(8, 8)) as single:
            with self.assertRaises(ValueError):
                single.temperature_grid("FIR")

    def test_temperature_node_weight_policy(self):
        with Preparation(GALAXIES, Grid(2, 2)) as p:
            tables = []
            try:
                # Zero stellar/FIR contributions isolate the CMB blend.
                for value in (0,0,0,0,2,10,0,0):
                    handle=ct.c_void_p()
                    check(p._lib, p._lib.cg_table_import(2,2,doubles([1,2]),
                          doubles([1,2]),doubles([value]*4),ct.byref(handle)))
                    tables.append(Table(p._lib,handle,{}))
                handle=ct.c_void_p()
                check(p._lib, p._lib.cg_combine_ic(doubles(GALAXIES[0].as_row()),
                      (ct.c_void_p*8)(*(t._handle for t in tables)),0.,0.,ct.byref(handle)))
                with Table(p._lib,handle,{}) as result:
                    self.assertEqual(result.snapshot().values,(10.,)*4)
                # The legacy upper endpoint has no next plane: explicit rejection.
                self.assertEqual(p._lib.cg_combine_ic(doubles(GALAXIES[0].as_row()),
                    (ct.c_void_p*8)(*(t._handle for t in tables)),1.,0.,ct.byref(handle)),1)
            finally:
                for table in tables:
                    table.close()

    def test_numpy_ownership_and_views(self):
        with Preparation(GALAXIES, Grid(8, 8)) as p:
            table = p.table("bs")
            self.assertEqual(table.values.shape, (8, 8))
            self.assertEqual(table.values.dtype, np.dtype("float64"))
            self.assertTrue(table.values.flags.c_contiguous)
            for buffer in (table._x, table._y, table._values):
                self.assertTrue(buffer.flags.owndata)
                self.assertFalse(buffer.flags.writeable)
            with self.assertRaises(ValueError):
                table.values[0, 0] = 1.
            with self.assertRaises(ValueError):
                table.values.flags.writeable = True
            retained = table.values
            expected = retained.copy()
        gc.collect()
        # Closing C descriptor must not free retained NumPy storage.
        np.testing.assert_array_equal(retained, expected)

    def test_cache_import_numpy_lifetime(self):
        with Preparation(GALAXIES, Grid(2, 2)) as p:
            x = np.array([1., 2.])
            y = np.array([1., 2.])
            z = np.array([2., 3., 3., 4.])
            table = Table.from_data(p._lib, TableData(x, y, z), {})
            x[:] = 99.
            del x, y, z
            gc.collect()
            with table:
                self.assertEqual(table.evaluate([1.5], [1.5]), [3.])
                self.assertTrue(table._values.flags.owndata)

    def test_native_borrow_is_zero_copy(self):
        with Preparation(GALAXIES, Grid(2, 2)) as p:
            x = np.array([1., 2.], dtype=np.float64)
            y = np.array([1., 2.], dtype=np.float64)
            z = np.array([2., 3., 3., 4.], dtype=np.float64)
            ptr = Table._pointer
            handle = ct.c_void_p()
            check(p._lib, p._lib.cg_table_borrow(2, 2, ptr(x), ptr(y), ptr(z),
                  ct.byref(handle)))
            try:
                out = doubles([0.])
                check(p._lib, p._lib.cg_table_eval(handle, 1, doubles([1.5]),
                      doubles([1.5]), out))
                self.assertEqual(out[0], 3.)
                # Test-only sequential mutation proves C borrowed the original
                # arrays. The public Python Table prevents such mutation.
                z[:] = 7.
                check(p._lib, p._lib.cg_table_eval(handle, 1, doubles([1.5]),
                      doubles([1.5]), out))
                self.assertEqual(out[0], 7.)
            finally:
                p._lib.cg_table_destroy(handle)
            np.testing.assert_array_equal(z, [7.] * 4)

if __name__ == "__main__":
    unittest.main()
