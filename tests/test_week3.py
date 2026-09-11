"""All-row comparisons against the independently built production driver."""
import ctypes as ct
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from congruents import Preparation, Grid, SolverGrid, load_catalogue, solve
from congruents.solver import _load, _TableInput

ROOT = Path(__file__).resolve().parents[1]
FILES = (
    ("transport", ("fcal","Dp","Dd","Dh","Q1","Q2","protons")),
    ("electrons", ("e1d","e2d","e1h","e2h")),
    ("emission", ("ic1d","ic2d","bs1d","bs2d","sy1d","sy2d","ic1h","ic2h",
                  "sy1h","sy2h","ff","tau_ff","pi","pi_fcal1","nu")),
)


class Week3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.folder = Path(cls.tmp.name)
        cls.catalogue = load_catalogue(ROOT/"input/cat_nt.txt")
        cls.p = Preparation(cls.catalogue, Grid(16,16))
        cls.addClassCleanup(cls.p.close)
        cls.grid = SolverGrid(16,8,16)
        cls.result = solve(cls.p, cls.grid, threads=4)

    def reference(self, full=False, medium=False):
        folder = self.folder/("medium" if medium else "full" if full else "legacy")
        data, output = folder/"data", folder/"output"
        if not output.exists():
            data.mkdir(parents=True)
            (output/"tau_loss").mkdir(parents=True)
            env = {**os.environ, "OMP_NUM_THREADS": "4"}
            commands = [
                [str(ROOT/"build"/("reference_precompute_medium" if medium else
                                  "reference_precompute_full" if full else "reference_precompute")),
                 str(ROOT/"input/cat_nt.txt"), str(data)],
                [str(ROOT/"build"/("reference_spectra_medium" if medium else "reference_spectra")),
                 str(ROOT/"input/cat_nt.txt"), str(data), str(output)],
            ]
            for i, command in enumerate(commands):
                with (folder/f"reference-{i}.log").open("w") as log:
                    subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                   stderr=subprocess.STDOUT, check=True, timeout=180)
        return output

    def compare(self, result, folder):
        for group, files in FILES:
            arrays = getattr(result, group)
            self.assertEqual(len(arrays), len(files))
            for (name, got), file in zip(arrays.items(), files):
                with self.subTest(component=name):
                    expected = np.loadtxt(folder/(file+".txt"), skiprows=1)
                    self.assertEqual(got.shape, expected.shape)
                    self.assertEqual(got.shape[0], 11)
                    self.assertTrue(np.isfinite(got).all())
                    self.assertTrue((got >= 0).all())
                    np.testing.assert_allclose(got, expected, rtol=2e-6, atol=1e-280)

    def test_all_components_all_galaxies_legacy_precision(self):
        self.compare(self.result, self.reference())

    def test_full_precision_separate_comparison(self):
        full = solve(self.p, self.grid, threads=4, legacy_table_precision=False)
        self.compare(full, self.reference(full=True))

    def test_denser_grid_all_components_all_galaxies(self):
        with Preparation(self.catalogue, Grid(32,32)) as preparation:
            result = solve(preparation, SolverGrid(32,16,32), threads=4)
        self.compare(result, self.reference(medium=True))

    def test_thread_count_equivalence(self):
        serial = solve(self.p, self.grid, threads=1)
        for group, _ in FILES:
            for name, values in getattr(self.result, group).items():
                np.testing.assert_array_equal(values, getattr(serial, group)[name])

    def test_save_and_readonly(self):
        path = self.folder/"result.npz"
        self.result.save(path)
        with np.load(path, allow_pickle=False) as data:
            np.testing.assert_array_equal(data["electrons__primary_disc"],
                                          self.result.electrons["primary_disc"])
            metadata = json.loads(str(data["metadata_json"]))
            self.assertTrue(metadata["legacy_table_precision"])
            self.assertEqual(metadata["solver_cells"], 16)
            self.assertEqual(len(metadata["catalogue"]), 11)
        with self.assertRaises(ValueError):
            self.result.electrons["primary_disc"][0,0] = 0.

    def test_validation_and_recovery(self):
        for arguments in ((3,8,16), (16,8,501), (16.,8,16)):
            with self.assertRaises((TypeError,ValueError)):
                SolverGrid(*arguments)
        with self.assertRaises(ValueError):
            solve(self.p, self.grid, threads=0)
        with self.assertRaises(FileNotFoundError):
            solve(self.p, self.grid, library=self.folder/"absent.dylib")
        # Invalid native buffers fail before entering any OpenMP worker.
        lib = _load()
        self.assertEqual(lib.cg_transport_batch(1,0,0,None,None,None,None,None,None), 1)
        self.assertEqual(lib.cg_solver_batch(1,0,0,0,0,
                         *([None]*7),*([None]*4),None,None,None), 1)
        # The coarse-grid halo solve produces non-finite integrands on this
        # reference catalogue. Fail safely, then prove a valid run still works.
        with Preparation(self.catalogue, Grid(8,8)) as p:
            with self.assertRaisesRegex(RuntimeError, "catalogue indices"):
                solve(p, SolverGrid(8,8,8), threads=4, legacy_table_precision=False)
        recovered = solve(self.p,self.grid,threads=4)
        np.testing.assert_array_equal(recovered.electrons["primary_disc"],
                                      self.result.electrons["primary_disc"])

    def test_private_gsl_and_serial_generation_not_exported(self):
        lib = _load()
        for name in ("gsl_set_error_handler", "init_do_2D_IC", "CRe_steadystate_solve"):
            with self.assertRaises(AttributeError):
                getattr(lib,name)

    def test_solver_equations_unchanged(self):
        source = (ROOT/"CRe_steadystate.h").read_text()
        source = source[:source.index("int CRe_steadystate_solve_number")]
        start = source.index("    int solve_system(")
        end = source.index("{", start)+1
        depth = 1
        while depth:
            depth += (source[end] == "{") - (source[end] == "}")
            end += 1
        expected = (source[:start]+source[end:]).replace("solve_system(", "cg_linear_solve(")
        expected = expected[expected.index("struct F_int_data"):].strip()
        adapter = (ROOT/"csrc/steady_state_native.h").read_text()
        actual = adapter[adapter.index("struct F_int_data"):adapter.rindex("#endif")].strip()
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
