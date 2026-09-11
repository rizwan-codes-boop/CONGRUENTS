"""Portable checks requiring no legacy GSL/cubature reference executable."""
import math
import unittest
from pathlib import Path
import numpy as np
from congruents import Preparation, Grid, load_catalogue
from congruents.preparation import Table, TableData
from congruents.quadrature import integrate
from congruents import serial

ROOT = Path(__file__).resolve().parents[1]


class SerialTests(unittest.TestCase):
    def test_quadrature_analytic(self):
        for f, low, high, expected in (
            (lambda x: x**4, 0., 1., .2),
            (math.sin, 0., math.pi, 2.),
            (lambda x: math.exp(-100*x), 0., 1., .01),
            (lambda x: 0., -1., 1., 0.),
        ):
            self.assertAlmostEqual(integrate(f, low, high), expected, places=12)
        with self.assertRaises(ValueError):
            integrate(math.sin, 1., 0.)
        with self.assertRaises(RuntimeError):
            integrate(lambda x: float("nan"), 0., 1.)

    def test_table_interpolation(self):
        with Table(TableData([1.,2.], [1.,2.], [2.,3.,3.,4.]), {}) as t:
            self.assertEqual(t.evaluate([1.,1.5,2.,3.], [1.,1.5,2.,1.]), [2.,3.,4.,0.])
            with self.assertRaises(ValueError):
                t.evaluate([1.], [1.,2.])
        with Table(TableData([1.,2.], [1.], [2.,4.]), {}) as t:
            self.assertEqual(t.evaluate([0.,1.5,3.]), [0.,3.,0.])
        for data in (TableData([2.,1.],[1.],[1.,2.]),
                     TableData([1.,2.],[1.],[-1.,2.]),
                     TableData([1.,2.],[1.],[float("nan"),2.])):
            with self.assertRaises(ValueError):
                Table(data, {})

    def test_ionisation_reference_values(self):
        expected = [-2.0513467838531972e-16, -3.7936953668288616e-16,
                    -4.955261088812638e-16]
        np.testing.assert_allclose(serial.ionisation([.001, 1., 100.], 1.), expected,
                                   rtol=1e-14, atol=0.)

    def test_galaxy_parallel_consistency(self):
        galaxies = load_catalogue(ROOT/"input/cat_nt.txt")
        with Preparation(galaxies, Grid(8,8), threads=1) as a, \
             Preparation(galaxies, Grid(8,8), threads=4) as b:
            self.assertTrue(a._context.openmp_enabled)
            self.assertEqual(a.properties, b.properties)
            self.assertEqual(a.config, b.config)
            self.assertEqual(a.temperature_bounds, b.temperature_bounds)
            fixture = np.loadtxt(ROOT/"tests/fixtures/gal_data.txt", skiprows=1)
            for i, p in enumerate(a.properties):
                np.testing.assert_allclose(list(p.values())[:9], fixture[i], rtol=6e-7, atol=0.)
                factor = 3. if math.log10(galaxies[i].sfr_msun_per_year /
                                           galaxies[i].stellar_mass_msun) > -10 else 1.5
                self.assertEqual(p["halo_magnetic_field_gauss"], p["magnetic_field_gauss"]/factor)


if __name__ == "__main__":
    unittest.main()
