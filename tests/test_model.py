from __future__ import annotations

import struct
import unittest

from elara_x_nrlmsis import parameters as p
from elara_x_nrlmsis.constants import (
    cintann, cmag, csfx, csfxmod, ctide, cut, dmissing, itb0, itex, itgb0,
    nmag, nspec, nsplNO, nsplO1, nut, specmass,
)
from elara_x_nrlmsis.model import msiscalc, _reset_cache_for_testing
from elara_x_nrlmsis.utilities import alt2gph


def _bits(x: float) -> bytes:
    return struct.pack('>d', float(x))


def _seed() -> None:
    p.initparmspace()
    for subset in (p.TN,p.N2,p.O2,p.O1,p.HE,p.H1,p.AR,p.N1,p.OA,p.NO):
        subset.beta.fill(0.0)
    p.swg[:] = [True] * len(p.swg)
    p.smod[:] = [False] * len(p.smod)
    p.N2Rflag = False

    for ix in range(itb0):
        p.TN.beta[0,ix] = 0.00445 - 0.000022 * float(ix)
    p.TN.beta[0,itb0]=445.0; p.TN.beta[0,itgb0]=20.0; p.TN.beta[0,itex]=1200.0
    p.TN.beta[1,itex]=0.75; p.TN.beta[cintann,itex]=0.50; p.TN.beta[ctide,itex]=0.25
    p.TN.beta[csfx,itex]=0.015; p.TN.beta[1,itb0]=0.050; p.TN.beta[cintann,itb0]=0.025
    p.TN.beta[csfxmod,itex]=0.002; p.TN.beta[csfxmod+1,itex]=-0.001; p.TN.beta[csfxmod+2,itex]=0.0005
    p.smod[itex]=True
    for j in range(cmag,cmag+nmag): p.TN.beta[j,itex]=0.0
    p.TN.beta[cmag,itex]=0.020; p.TN.beta[cmag+1,itex]=0.0002; p.TN.beta[cmag+8,itex]=0.00001
    p.TN.beta[cmag+27,itex]=0.002; p.TN.beta[cmag+28,itex]=0.00001
    for j in range(cut,cut+nut): p.TN.beta[j,itex]=0.0
    p.TN.beta[cut,itex]=0.002; p.TN.beta[cut+1,itex]=-0.001

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

    p.specflag[:] = [True]*(nspec-1)
    p.massflag[:] = [True]*(nspec-1)
    p.masswgt[:] = list(specmass)
    p.masswgt[0]=0.0; p.masswgt[9]=0.0
    p.zaltflag=False
    p.initflag=True
    _reset_cache_for_testing()


class TestMainModel(unittest.TestCase):
    def setUp(self):
        _seed()
        self.ap = [4.0]*7

    def test_ap_requires_seven_values(self):
        with self.assertRaises(ValueError):
            msiscalc(100,43200,100,45,10,150,150,[4.0]*6)

    def test_default_return_shape(self):
        tn,dn = msiscalc(100,43200,100,45,10,150,150,self.ap)
        self.assertIsInstance(tn,float)
        self.assertEqual(len(dn),10)

    def test_optional_tex_return_shape(self):
        tn,dn,tex = msiscalc(100,43200,100,45,10,150,150,self.ap,return_tex=True)
        self.assertIsInstance(tex,float)
        self.assertGreater(tex,tn)
        self.assertEqual(len(dn),10)

    def test_identical_repeat_is_bitwise_equal(self):
        a = msiscalc(172.25,30000,80,20,-30,145,155,self.ap,return_tex=True)
        b = msiscalc(172.25,30000,80,20,-30,145,155,self.ap,return_tex=True)
        self.assertEqual(_bits(a[0]),_bits(b[0])); self.assertEqual(tuple(map(_bits,a[1])),tuple(map(_bits,b[1]))); self.assertEqual(_bits(a[2]),_bits(b[2]))

    def test_altitude_only_repeat_is_bitwise_equal(self):
        msiscalc(172.25,30000,80,20,-30,145,155,self.ap,return_tex=True)
        a = msiscalc(172.25,30000,90,20,-30,145,155,self.ap,return_tex=True)
        b = msiscalc(172.25,30000,90,20,-30,145,155,self.ap,return_tex=True)
        self.assertEqual(tuple(map(_bits,(a[0],*a[1],a[2]))),tuple(map(_bits,(b[0],*b[1],b[2]))))

    def test_return_to_baseline_is_bitwise_equal(self):
        a = msiscalc(172.25,30000,80,20,-30,145,155,self.ap,return_tex=True)
        msiscalc(173.25,36000,90,25,-20,150,180,[30,20,18,16,14,12,10],return_tex=True)
        b = msiscalc(172.25,30000,80,20,-30,145,155,self.ap,return_tex=True)
        self.assertEqual(tuple(map(_bits,(a[0],*a[1],a[2]))),tuple(map(_bits,(b[0],*b[1],b[2]))))

    def test_geodetic_and_equivalent_geopotential_are_bitwise_equal(self):
        p.zaltflag=True
        a=msiscalc(220,54000,100,45,25,160,170,self.ap,return_tex=True)
        zeta=alt2gph(45.0,100.0)
        p.zaltflag=False
        b=msiscalc(220,54000,zeta,45,25,160,170,self.ap,return_tex=True)
        self.assertEqual(tuple(map(_bits,(a[0],*a[1],a[2]))),tuple(map(_bits,(b[0],*b[1],b[2]))))

    def test_optional_tex_does_not_change_common_outputs(self):
        a=msiscalc(250.5,70000,140,-33.25,179,180,220,self.ap)
        b=msiscalc(250.5,70000,140,-33.25,179,180,220,self.ap,return_tex=True)
        self.assertEqual(_bits(a[0]),_bits(b[0])); self.assertEqual(tuple(map(_bits,a[1])),tuple(map(_bits,b[1])))

    def test_disabled_species_return_missing(self):
        p.specflag[4]=False; p.specflag[5]=False; p.massflag[4]=False; p.massflag[5]=False; p.masswgt[4]=0.; p.masswgt[5]=0.
        _,dn=msiscalc(301,40000,100,10,40,150,160,self.ap)
        self.assertEqual(dn[4],dmissing); self.assertEqual(dn[5],dmissing)

    def test_mass_density_gate_returns_missing(self):
        p.specflag[0]=False; p.massflag[:]=[False]*(nspec-1); p.masswgt[:]=[0.0]*(nspec-1)
        _,dn=msiscalc(302,40000,100,10,40,150,160,self.ap)
        self.assertEqual(dn[0],dmissing)

    def test_zetaF_equality_is_finite(self):
        tn,dn=msiscalc(100,43200,70.0,45,10,150,150,self.ap)
        self.assertTrue(all(x==x for x in (tn,*dn)))

    def test_zetaB_equality_is_finite(self):
        tn,dn=msiscalc(100,43200,122.5,45,10,150,150,self.ap)
        self.assertTrue(all(x==x for x in (tn,*dn)))

    def test_self_initialization_delegates_to_parameters(self):
        original=p.msisinit
        called=[]
        def fake_msisinit(*args,**kwargs):
            called.append(True); p.initflag=True
        try:
            p.initflag=False; p.msisinit=fake_msisinit
            msiscalc(100,43200,100,45,10,150,150,self.ap)
            self.assertEqual(called,[True])
        finally:
            p.msisinit=original; p.initflag=True


if __name__ == '__main__':
    unittest.main()
