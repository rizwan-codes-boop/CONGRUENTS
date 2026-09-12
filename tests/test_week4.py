"""Initial observer-output stage: geometry and explicit-tau transformation."""
from pathlib import Path
import subprocess
import unittest
import tempfile
import os
import sys
import numpy as np
from congruents import load_catalogue
from congruents.observer import luminosity_distance_mpc, distance_factor_cm2, observer_sed
from congruents import Preparation, Grid, SolverGrid
from congruents.pipeline import run, COMPONENT_FILES
from congruents.attenuation import EBLTable
from compare_outputs import compare

ROOT = Path(__file__).resolve().parents[1]


class ObserverTests(unittest.TestCase):
    def test_ebl_extrapolation_against_gsl(self):
        energies = [1e-16,1.,12.,1e3,1e6,1e8]
        for name in ("Franceschini","Dominguez","Gilmore"):
            table = EBLTable(ROOT/"input"/f"tau_Eg_z_{name}.txt")
            for z in (0.,.0001,.1,4.):
                expected = np.fromstring(subprocess.check_output([str(ROOT/"build/direct_observer"),
                    str(table.path),str(z),*map(str,energies)],text=True),sep=" ")
                np.testing.assert_allclose(table.optical_depth(energies,[z])[0],expected,rtol=1e-12,atol=1e-12)
    def test_end_to_end_all_catalogue_outputs(self):
        self._compare_run(16,8,"")

    def test_denser_end_to_end_all_catalogue_outputs(self):
        self._compare_run(32,16,"_medium")

    def _compare_run(self,cells,photons,suffix):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            data, output = folder/"data",folder/"reference"
            data.mkdir()
            (output/"tau_loss").mkdir(parents=True)
            for executable in ("reference_precompute"+suffix,"reference_observer"+suffix):
                args = [str(ROOT/"build"/executable),str(ROOT/"input/cat_nt.txt"),str(data)]
                if executable.startswith("reference_observer"):
                    args.append(str(output))
                with (folder/(executable+".log")).open("w") as log:
                    subprocess.run(args,cwd=ROOT,env={**os.environ,"OMP_NUM_THREADS":"4"},
                                   stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
            with Preparation(load_catalogue(ROOT/"input/cat_nt.txt"),Grid(cells,cells)) as p:
                result = run(p,SolverGrid(cells,photons,cells),threads=4)
                if cells==16:
                    single = run(p,SolverGrid(cells,photons,cells),threads=1)
                    for name in result.components:
                        np.testing.assert_array_equal(result.components[name],single.components[name])
                    for name in result.source.diagnostics:
                        np.testing.assert_array_equal(result.source.diagnostics[name],single.source.diagnostics[name])
            for name,file in COMPONENT_FILES.items():
                with self.subTest(component=name):
                    np.testing.assert_allclose(result.components[name],
                        np.loadtxt(output/file,skiprows=1),rtol=3e-6,atol=1e-280)
            np.testing.assert_allclose(result.tau_internal,np.loadtxt(output/"tau_gg.txt",skiprows=1),rtol=2e-6,atol=1e-280)
            np.testing.assert_allclose(result.gamma_luminosity,np.loadtxt(output/"L_gamma.txt",skiprows=1),rtol=3e-6)
            for name,array in result.source.diagnostics.items():
                file = output/"tau_loss"/(name+".txt") if name.startswith("tau_loss") else output/(name+".txt")
                with self.subTest(diagnostic=name):
                    expected = np.loadtxt(file,skiprows=0 if name.startswith("CR_specs") else 1)
                    np.testing.assert_allclose(array,expected,rtol=3e-6,atol=1e-280)
            result.export(folder/"python")
            report = compare(output,folder/"python")
            self.assertEqual({name:r for name,r in report.items() if not r["passed"]},{})
            with self.assertRaises(FileExistsError):
                result.export(folder/"python")

    def test_cli_refuses_existing_directory_before_work(self):
        with tempfile.TemporaryDirectory() as temp:
            proc = subprocess.run([sys.executable,"-m","congruents",str(ROOT/"input/cat_nt.txt"),temp],
                cwd=ROOT,env={**os.environ,"PYTHONPATH":str(ROOT/"src")},capture_output=True,text=True)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn("Output already exists",proc.stderr)

    def test_ebl_nodes_and_invalid_file(self):
        table = EBLTable(ROOT/"input/tau_Eg_z_Franceschini.txt")
        for i,z in enumerate(table.z):
            got = table.optical_depth(table.energy/1e9*(1+z),[z])[0]
            np.testing.assert_allclose(got,table.values[i],rtol=1e-12,atol=1e-12)
        self.assertTrue(np.all(table.optical_depth([1e-16,1e8],[0.,.001,4.])>=0))
        with tempfile.NamedTemporaryFile(mode="w") as f:
            f.write("2 2 0 1");f.flush()
            with self.assertRaises(ValueError):
                EBLTable(f.name)

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
