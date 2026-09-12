"""Initial observer-output stage: geometry and explicit-tau transformation."""
from pathlib import Path
import subprocess
import unittest
import numpy as np
from congruents import load_catalogue
from congruents.observer import luminosity_distance_mpc, distance_factor_cm2, observer_sed

ROOT = Path(__file__).resolve().parents[1]


class ObserverTests(unittest.TestCase):
    def test_reference_geometry_and_transform_all_galaxies(self):
        rows = np.array([g.as_row() for g in load_catalogue(ROOT/"input/cat_nt.txt")])
        z = np.r_[rows[:,0], .1, 1., 3.]
        expected = np.array([np.fromstring(subprocess.check_output(
            [str(ROOT/"build/direct_observer"),str(v)],text=True),sep=" ") for v in z])
        np.testing.assert_allclose(luminosity_distance_mpc(z),expected[:,0],rtol=1e-12)
        np.testing.assert_allclose(distance_factor_cm2(z),expected[:,1],rtol=1e-12)
        shape = (len(z),4)
        source = np.broadcast_to([3.,7.,11.,19.],shape)
        ti = np.broadcast_to(np.arange(4)*.1,shape)
        te = np.broadcast_to(np.arange(4)*.2,shape)
        result = observer_sed([1.,2.,4.,8.],source,z,tau_internal=ti,tau_ebl=te)
        np.testing.assert_allclose(result,expected[:,2:],rtol=1e-12,atol=0.)
        np.testing.assert_array_equal(result[:,-1],0.)
        self.assertFalse(result.flags.writeable)
        np.testing.assert_array_equal(source[0],[3.,7.,11.,19.])

    def test_required_attenuation_and_validation(self):
        for z in ([0.],[-1.],[np.nan],[],.1):
            with self.assertRaises(ValueError):
                luminosity_distance_mpc(z)
        with self.assertRaises(TypeError):
            observer_sed([1,2],[[1,1]],[.1])
        for e, source, tau in (([2,1],[[1,1]],[[0,0]]),
                              ([1,2],[[1,-1]],[[0,0]]),
                              ([1,2],[[1,1]],[[0,-1]]),
                              ([1,2],[[1,1]],[0,0])):
            with self.assertRaises(ValueError):
                observer_sed(e,source,[.1],tau_internal=tau,tau_ebl=[[0,0]])

    def test_explicit_zero_attenuation_diagnostic(self):
        got = observer_sed([1,2,4],[[1,2,4]],[1.],
                           tau_internal=[[0,0,0]],tau_ebl=[[0,0,0]])
        np.testing.assert_allclose(got[0],np.array([2,16,0])*distance_factor_cm2([1.])[0])
