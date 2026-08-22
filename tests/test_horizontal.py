from __future__ import annotations

import math
import struct
import unittest

from elara_x_nrlmsis import horizontal as h
from elara_x_nrlmsis import parameters as p
from elara_x_nrlmsis.constants import (
    cextra,
    cmag,
    cspw,
    csfx,
    csfxmod,
    ctide,
    cut,
    maxnbf,
    mbf,
    nmag,
    nut,
)


def f32(value: float) -> float:
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


class NativeHorizontalContractTests(unittest.TestCase):
    def setUp(self):
        p.swg[:] = [True] * maxnbf
        p.zsfx[:] = [False] * (mbf + 1)
        p.tsfx[:] = [False] * (mbf + 1)
        p.psfx[:] = [False] * (mbf + 1)
        h.lastlat = f32(-999.9)
        h.lastdoy = f32(-999.9)
        h.lastlst = f32(-999.9)
        h.lastlon = f32(-999.9)
        for array in (h.cdoy, h.sdoy, h.clon, h.slon, h.clst, h.slst):
            for i in range(len(array)):
                array[i] = 0.0
        for n in range(len(h.plg)):
            for m in range(len(h.plg[n])):
                h.plg[n][m] = 0.0

    def test_default_real_cache_sentinels_are_promoted_binary32(self):
        self.assertEqual(h.lastlat, -999.9000244140625)
        self.assertEqual(h.lastdoy, -999.9000244140625)
        self.assertEqual(h.lastlst, -999.9000244140625)
        self.assertEqual(h.lastlon, -999.9000244140625)

    def test_solzen_coefficients_preserve_default_real_rounding(self):
        expected = tuple(f32(x) for x in (0.017203534, 0.034407068, 0.051610602, 0.068814136, 0.103221204))
        self.assertEqual(h._SOLZEN_P, expected)

    def test_solzen_authoritative_longitude_overwrite_is_preserved(self):
        a = h.solzen(80.0, 12.0, 0.0, 0.0)
        b = h.solzen(80.0, 12.0, 0.0, 173.25)
        self.assertEqual(a, b)

    def test_globe_fortran_mod_semantics_can_produce_negative_lst(self):
        h.globe(100.0, -90000.0, 15.0, -180.0, 150.0, 150.0, [4.0] * 7)
        self.assertEqual(h.lastlst, -13.0)

    def test_globe_preserves_basis_layout_gap_at_cmag_plus_7(self):
        bf = h.globe(100.0, 43200.0, 45.0, 0.0, 150.0, 150.0, [4.0] * 7)
        self.assertEqual(len(bf), maxnbf)
        self.assertEqual(bf[cmag + 7], 0.0)

    def test_globe_final_switch_mask_only_zeroes_linear_range(self):
        p.swg[:] = [False] * maxnbf
        bf = h.globe(100.0, 43200.0, 45.0, 0.0, 175.0, 190.0, [4.0] * 7)
        self.assertTrue(all(x == 0.0 for x in bf[: mbf + 1]))
        self.assertEqual(bf[csfxmod], 25.0)

    def test_globe_in_place_output_contract(self):
        bf = [99.0] * maxnbf
        returned = h.globe(100.0, 43200.0, 45.0, 0.0, 150.0, 150.0, [4.0] * 7, bf)
        self.assertIs(returned, bf)
        self.assertEqual(bf[0], 1.0)

    def test_sfluxmod_overlap_uses_zonal_then_cycle_precedence(self):
        subset = p.BasisSubset(name="TEST", bl=0, nl=0, beta=p.FortranMatrix(maxnbf, 0, 0, 0.0))
        gf = [0.0] * maxnbf
        subset.beta[csfxmod, 0] = 0.25
        gf[csfxmod] = 0.75
        subset.beta[9, 0] = 0.5
        gf[9] = 0.25
        p.zsfx[9] = p.tsfx[9] = p.psfx[9] = True
        value = h.sfluxmod(0, gf, subset, 1.0)
        expected = 0.5 * 0.25 * (0.25 * 0.75)
        self.assertEqual(value, expected)

    def test_geomag_returns_zero_when_both_master_switches_off(self):
        p.swg[cmag] = False
        p.swg[cmag + 1] = False
        self.assertEqual(h.geomag([0.0] * nmag, [0.0] * 13, [[0.0, 0.0] for _ in range(7)]), 0.0)

    def test_geomag_daily_zero_k00s_branch(self):
        p.swg[cmag] = True
        p.swg[cmag + 1] = True
        p0 = [0.0] * nmag
        p0[0] = 1.0
        p0[1] = 0.0
        self.assertEqual(h.geomag(p0, [0.0] * 13, [[0.0, 0.0] for _ in range(7)]), 0.0)

    def test_geomag_history_zero_beta_branch(self):
        p.swg[cmag] = True
        p.swg[cmag + 1] = False
        p0 = [0.0] * nmag
        p0[28] = 0.0
        self.assertEqual(h.geomag(p0, [0.0] * 13, [[0.0, 0.0] for _ in range(7)]), 0.0)

    def test_utdep_phase_parameters_are_not_masked(self):
        p0 = [0.0] * nut
        p0[0], p0[1], p0[2] = 0.25, -0.5, 0.75
        bf = [-0.5, -0.375, -0.25, -0.125, 0.25, 0.125, 0.25, 0.375, 0.5]
        p.swg[cut : cut + nut] = [False] * nut
        value = h.utdep(p0, bf)
        self.assertEqual(value, 0.0)
        # Turning on one amplitude term must use the still-present phase p(0).
        p0[6] = 1.0
        p.swg[cut + 6] = True
        expected = math.cos(bf[0] - p0[0]) * (p0[6] * bf[4])
        self.assertEqual(h.utdep(p0, bf), expected)

    def test_input_shape_guards_follow_fortran_array_contracts(self):
        with self.assertRaises(ValueError):
            h.globe(1, 0, 0, 0, 150, 150, [4.0] * 6)
        with self.assertRaises(ValueError):
            h.geomag([0.0] * (nmag - 1), [0.0] * 13, [[0.0, 0.0] for _ in range(7)])
        with self.assertRaises(ValueError):
            h.utdep([0.0] * nut, [0.0] * 8)


if __name__ == "__main__":
    unittest.main()
