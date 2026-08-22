from __future__ import annotations

import struct
import unittest

from elara_x_nrlmsis import legacy_interface as li
from elara_x_nrlmsis import parameters as p
from elara_x_nrlmsis.constants import (
    cintann, cmag, csfx, csfxmod, ctide, cut, dmissing, itb0, itex, itgb0,
    nmag, nspec, nsplNO, nsplO1, nut, specmass,
)
from elara_x_nrlmsis.model import _reset_cache_for_testing


def _bits32(x: float) -> bytes:
    return struct.pack(">f", float(x))


def _f32(x: float) -> float:
    return struct.unpack("<f", struct.pack("<f", float(x)))[0]


def _seed() -> None:
    p.initparmspace()
    for subset in (p.TN, p.N2, p.O2, p.O1, p.HE, p.H1, p.AR, p.N1, p.OA, p.NO):
        subset.beta.fill(0.0)
    p.swg[:] = [True] * len(p.swg)
    p.smod[:] = [False] * len(p.smod)
    p.N2Rflag = False

    for ix in range(itb0):
        p.TN.beta[0, ix] = 0.00445 - 0.000022 * float(ix)
    p.TN.beta[0, itb0] = 445.0
    p.TN.beta[0, itgb0] = 20.0
    p.TN.beta[0, itex] = 1200.0
    p.TN.beta[1, itex] = 0.75
    p.TN.beta[cintann, itex] = 0.50
    p.TN.beta[ctide, itex] = 0.25
    p.TN.beta[csfx, itex] = 0.015
    p.TN.beta[1, itb0] = 0.050
    p.TN.beta[cintann, itb0] = 0.025
    p.TN.beta[csfxmod, itex] = 0.002
    p.TN.beta[csfxmod + 1, itex] = -0.001
    p.TN.beta[csfxmod + 2, itex] = 0.0005
    p.smod[itex] = True
    for j in range(cmag, cmag + nmag):
        p.TN.beta[j, itex] = 0.0
    p.TN.beta[cmag, itex] = 0.020
    p.TN.beta[cmag + 1, itex] = 0.0002
    p.TN.beta[cmag + 8, itex] = 0.00001
    p.TN.beta[cmag + 27, itex] = 0.002
    p.TN.beta[cmag + 28, itex] = 0.00001
    for j in range(cut, cut + nut):
        p.TN.beta[j, itex] = 0.0
    p.TN.beta[cut, itex] = 0.002
    p.TN.beta[cut + 1, itex] = -0.001

    p.N2.beta[0,1]=100.; p.N2.beta[0,2]=8.; p.N2.beta[0,3]=15.; p.N2.beta[0,7]=.04; p.N2.beta[0,8]=105.; p.N2.beta[0,9]=24.
    p.O2.beta[0,1]=101.; p.O2.beta[0,2]=8.5; p.O2.beta[0,3]=15.5; p.O2.beta[0,7]=.035; p.O2.beta[0,8]=104.; p.O2.beta[0,9]=23.
    p.O1.beta[0,0]=36.; p.O1.beta[0,1]=100.; p.O1.beta[0,2]=8.; p.O1.beta[0,3]=15.; p.O1.beta[0,4]=.15; p.O1.beta[0,5]=92.; p.O1.beta[0,6]=20.; p.O1.beta[0,7]=.045; p.O1.beta[0,8]=105.; p.O1.beta[0,9]=25.
    for j in range(nsplO1): p.O1.beta[0,j+10]=34.8-.18*float(j)
    p.HE.beta[0,1]=102.; p.HE.beta[0,2]=9.; p.HE.beta[0,3]=16.; p.HE.beta[0,7]=.03; p.HE.beta[0,8]=108.; p.HE.beta[0,9]=28.
    p.H1.beta[0,0]=29.5; p.H1.beta[0,1]=100.; p.H1.beta[0,2]=8.; p.H1.beta[0,3]=16.; p.H1.beta[0,4]=.12; p.H1.beta[0,5]=95.; p.H1.beta[0,6]=22.; p.H1.beta[0,7]=.04; p.H1.beta[0,8]=108.; p.H1.beta[0,9]=28.
    p.AR.beta[0,1]=99.; p.AR.beta[0,2]=8.; p.AR.beta[0,3]=14.; p.AR.beta[0,7]=.025; p.AR.beta[0,8]=102.; p.AR.beta[0,9]=22.
    p.N1.beta[0,0]=30.; p.N1.beta[0,1]=101.; p.N1.beta[0,2]=8.; p.N1.beta[0,3]=16.; p.N1.beta[0,4]=.11; p.N1.beta[0,5]=96.; p.N1.beta[0,6]=22.; p.N1.beta[0,7]=.035; p.N1.beta[0,8]=107.; p.N1.beta[0,9]=27.
    p.OA.beta[0,0]=25.; p.OA.beta[0,4]=.08; p.OA.beta[0,5]=180.; p.OA.beta[0,6]=120.
    p.NO.beta[0,0]=31.; p.NO.beta[0,1]=101.; p.NO.beta[0,2]=8.; p.NO.beta[0,3]=16.; p.NO.beta[0,4]=.10; p.NO.beta[0,5]=96.; p.NO.beta[0,6]=21.; p.NO.beta[0,7]=.04; p.NO.beta[0,8]=108.; p.NO.beta[0,9]=28.
    for j in range(nsplNO): p.NO.beta[0,j+10]=30.8-.14*float(j)
    p.O1.beta[1,0]=.002; p.N1.beta[cintann,0]=.002; p.HE.beta[csfx,7]=.0005; p.NO.beta[1,0]=.002

    p.specflag[:] = [True] * (nspec - 1)
    p.massflag[:] = [True] * (nspec - 1)
    p.masswgt[:] = list(specmass)
    p.masswgt[0] = 0.0
    p.masswgt[9] = 0.0
    p.zaltflag = False
    p.initflag = True
    _reset_cache_for_testing()


class TestLegacyGTD8D(unittest.TestCase):
    def setUp(self):
        _seed()
        self.ap = [4.0] * 7

    def test_ap_requires_seven_values(self):
        with self.assertRaises(ValueError):
            li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, [4.0] * 6, 0)

    def test_fortran_mod_positive_and_negative(self):
        self.assertEqual(li._fortran_mod_int(99123, 1000), 123)
        self.assertEqual(li._fortran_mod_int(-10123, 1000), -123)
        self.assertEqual(li._fortran_mod_int(-123, 1000), -123)

    def test_return_shape_and_binary32_boundary(self):
        d, t = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 0)
        self.assertEqual(len(d), 10)
        self.assertEqual(len(t), 2)
        self.assertTrue(all(_f32(x) == x for x in (*d, *t)))

    def test_stl_is_ignored_bitwise(self):
        a = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 0)
        b = li.gtd8d(24172, 30000, 80, 20, -30, 23.75, 145, 155, self.ap, 0)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_mass_is_ignored_bitwise(self):
        a = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 0)
        b = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 999)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_year_digits_are_ignored(self):
        a = li.gtd8d(99123, 43200, 100, 45, 10, 5, 150, 150, self.ap, 48)
        b = li.gtd8d(123, 43200, 100, 45, 10, 5, 150, 150, self.ap, 48)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_negative_iyd_uses_fortran_mod_semantics(self):
        a = li.gtd8d(-10123, 43200, 100, 45, 10, 5, 150, 150, self.ap, 48)
        b = li.gtd8d(-123, 43200, 100, 45, 10, 5, 150, 150, self.ap, 48)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_float32_input_collapse_is_exact(self):
        a = li.gtd8d(24123, 43200.125, 100.000001, 45.123456, -123.98765, 6.75, 150.125, 151.5, self.ap, 48)
        b = li.gtd8d(24123, 43200.125, 100.000002, 45.123456, -123.98765, 6.75, 150.125, 151.5, self.ap, 48)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_species_reorder_and_unit_conversion(self):
        original = li.msiscalc
        try:
            native = (1.0e-9, 2.0e10, 3.0e10, 4.0e10, 5.0e10, 6.0e10, 7.0e10, 8.0e10, 9.0e10, 1.0e11)
            li.msiscalc = lambda *args, **kwargs: (250.25, native, 1000.5)
            d, t = li.gtd8d(24172, 1, 2, 3, 4, 5, 6, 7, self.ap, 8)
            expected = tuple(map(_f32, (5.0e4, 4.0e4, 2.0e4, 3.0e4, 7.0e4, 1.0e-12, 6.0e4, 8.0e4, 9.0e4, 1.0e5)))
            self.assertEqual(tuple(map(_bits32, d)), tuple(map(_bits32, expected)))
            self.assertEqual(tuple(map(_bits32, t)), tuple(map(_bits32, (_f32(1000.5), _f32(250.25)))))
        finally:
            li.msiscalc = original

    def test_missing_density_is_not_scaled(self):
        original = li.msiscalc
        try:
            native = (dmissing, dmissing, 3.0e10, 4.0e10, dmissing, 6.0e10, 7.0e10, 8.0e10, 9.0e10, 1.0e11)
            li.msiscalc = lambda *args, **kwargs: (250.0, native, 1000.0)
            d, _ = li.gtd8d(24172, 1, 2, 3, 4, 5, 6, 7, self.ap, 8)
            missing32 = _f32(dmissing)
            self.assertEqual(_bits32(d[0]), _bits32(missing32))
            self.assertEqual(_bits32(d[2]), _bits32(missing32))
            self.assertEqual(_bits32(d[5]), _bits32(missing32))
        finally:
            li.msiscalc = original

    def test_repeated_call_is_bitwise_identical(self):
        a = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 0)
        b = li.gtd8d(24172, 30000, 80, 20, -30, 0, 145, 155, self.ap, 0)
        self.assertEqual(tuple(map(_bits32, (*a[0], *a[1]))), tuple(map(_bits32, (*b[0], *b[1]))))

    def test_disabled_he_and_h_map_to_legacy_missing_slots(self):
        p.specflag[4] = False
        p.specflag[5] = False
        p.massflag[4] = False
        p.massflag[5] = False
        p.masswgt[4] = 0.0
        p.masswgt[5] = 0.0
        d, _ = li.gtd8d(24301, 40000, 100, 10, 40, 0, 150, 160, self.ap, 0)
        missing32 = _f32(dmissing)
        self.assertEqual(_bits32(d[0]), _bits32(missing32))
        self.assertEqual(_bits32(d[6]), _bits32(missing32))

    def test_float32_reachable_thresholds_are_finite(self):
        for altitude in (69.99999, 70.0, 70.00001, 122.49999, 122.5, 122.50001):
            d, t = li.gtd8d(24100, 43200, altitude, 45, 10, 0, 150, 150, self.ap, 0)
            self.assertTrue(all(x == x for x in (*d, *t)))


if __name__ == "__main__":
    unittest.main()
