from __future__ import annotations

import math
from pathlib import Path
import unittest

from elara_x_nrlmsis import constants as c
from elara_x_nrlmsis import parameters as p
from elara_x_nrlmsis import temperature as t
from elara_x_nrlmsis.utilities import bspline


class NativeTemperatureContractTests(unittest.TestCase):
    def setUp(self):
        p.initparmspace()
        p.swg[:] = [True] * len(p.swg)
        p.smod[:] = [False] * len(p.smod)

    def tearDown(self):
        p.swg[:] = [True] * len(p.swg)
        p.smod[:] = [False] * len(p.smod)

    def _baseline(self):
        gf = [0.0] * c.maxnbf
        gf[0] = 1.0
        p.TN.beta.fill(0.0)
        basecf = 0.00440
        slope = -0.000025
        for ix in range(0, c.itb0):
            p.TN.beta[0, ix] = basecf + slope * float(ix)
        p.TN.beta[0, c.itb0] = 440.0
        p.TN.beta[0, c.itgb0] = 21.0
        p.TN.beta[0, c.itex] = 1150.0
        return gf

    def _synthetic_tpro(self):
        s = t.TnParm()
        for j in range(c.nl + 1):
            s.cf[j] = 0.0025 + 0.00005 * float(j)
        s.tex = 1200.0
        s.tb0 = 420.0
        s.sigma = 0.04
        return s

    def test_tnparm_array_bounds(self):
        s = t.TnParm()
        self.assertEqual(len(s.cf), c.nl + 1)
        self.assertEqual(len(s.beta), c.nl + 1)
        self.assertEqual(len(s.gamma), c.nl + 1)

    def test_tfnparm_baseline_primary_outputs(self):
        s = t.tfnparm(self._baseline())
        self.assertEqual(s.tex, 1150.0)
        self.assertEqual(s.tgb0, 21.0)
        self.assertEqual(s.tb0, 440.0)

    def test_tfnparm_unconstrained_cf_order(self):
        s = t.tfnparm(self._baseline())
        for ix in range(0, c.itb0):
            self.assertEqual(s.cf[ix], 0.00440 - 0.000025 * float(ix))

    def test_tfnparm_c2_matmul_orientation(self):
        s = t.tfnparm(self._baseline())
        bc0 = 1.0 / s.tb0
        bc1 = -s.tgb0 / (s.tb0 * s.tb0)
        bc2 = -bc1 * (s.sigma + 2.0 * s.tgb0 / s.tb0)
        bc = (bc0, bc1, bc2)
        expected = [
            sum(bc[i] * c.c2tn[i][j] for i in range(3))
            for j in range(3)
        ]
        for j in range(3):
            self.assertAlmostEqual(s.cf[c.itb0 + j], expected[j], places=15)

    def test_beta_recurrence_is_sequential(self):
        s = t.tfnparm(self._baseline())
        running = s.cf[0] * c.wbeta[0]
        self.assertAlmostEqual(s.beta[0], running, places=15)
        for ix in range(1, c.nl + 1):
            running = running + s.cf[ix] * c.wbeta[ix]
            self.assertAlmostEqual(s.beta[ix], running, places=15)

    def test_gamma_recurrence_is_sequential(self):
        s = t.tfnparm(self._baseline())
        running = s.beta[0] * c.wgamma[0]
        self.assertAlmostEqual(s.gamma[0], running, places=15)
        for ix in range(1, c.nl + 1):
            running = running + s.beta[ix] * c.wgamma[ix]
            self.assertAlmostEqual(s.gamma[ix], running, places=15)

    def test_tfnparm_integration_scalars_are_finite(self):
        s = t.tfnparm(self._baseline())
        values = (
            s.tzetaF, s.tzetaA, s.dlntdzA, s.lndtotF, s.sigma, s.sigmasq,
            s.b, s.cVs, s.cVb, s.cWs, s.cWb, s.VzetaF, s.VzetaA,
            s.WzetaA, s.Vzeta0,
        )
        self.assertTrue(all(math.isfinite(x) for x in values))

    def test_tfnx_iz0_preserves_negative_logical_weights(self):
        s = self._synthetic_tpro()
        w = [0.10, 0.20, 0.30, 0.40]
        expected = 1.0 / (s.cf[0] * 0.40)
        self.assertEqual(t.tfnx(0.0, 0, w, s), expected)

    def test_tfnx_iz2_preserves_logical_slice(self):
        s = self._synthetic_tpro()
        w = [0.10, 0.20, 0.30, 0.40]
        expected = 1.0 / (
            s.cf[0] * 0.20 + s.cf[1] * 0.30 + s.cf[2] * 0.40
        )
        self.assertEqual(t.tfnx(20.0, 2, w, s), expected)

    def test_tfnx_mapping_weights_are_supported_by_logical_index(self):
        s = self._synthetic_tpro()
        w = {-3: 0.10, -2: 0.20, -1: 0.30, 0: 0.40}
        expected = 1.0 / sum(
            s.cf[i] * w[-3 + i] for i in range(4)
        )
        self.assertEqual(t.tfnx(30.0, 3, w, s), expected)

    def test_tfnx_exact_zetaB_uses_bates_branch(self):
        s = self._synthetic_tpro()
        self.assertEqual(t.tfnx(c.zetaB, 0, [0.0] * 4, s), s.tb0)

    def test_tfnx_physical_spline_to_bates_boundary(self):
        s = t.tfnparm(self._baseline())
        z = 122.499999
        spl, iz = bspline(z, c.nodesTN, c.nd, 4, p.etaTN)
        w = [spl[-3, 4], spl[-2, 4], spl[-1, 4], spl[0, 4]]
        below = t.tfnx(z, iz, w, s)
        at = t.tfnx(c.zetaB, 0, [0.0] * 4, s)
        self.assertLess(abs(below - at), 1.0e-3)

    def test_temperature_module_has_no_pymsis_or_fortran_runtime_dependency(self):
        text = Path(t.__file__).read_text(encoding="utf-8").lower()
        self.assertNotIn("import pymsis", text)
        self.assertNotIn("from pymsis", text)
        self.assertNotIn("subprocess", text)
        self.assertNotIn("ctypes", text)


if __name__ == "__main__":
    unittest.main()
