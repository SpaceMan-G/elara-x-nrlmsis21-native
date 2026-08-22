from __future__ import annotations

import math
import unittest

from elara_x_nrlmsis import constants as c
from elara_x_nrlmsis import density as d
from elara_x_nrlmsis import parameters as p
from elara_x_nrlmsis import temperature as t
from elara_x_nrlmsis import utilities as u


def _seed():
    p.initparmspace()
    p.swg[:] = [True] * len(p.swg)
    p.smod[:] = [False] * len(p.smod)
    gf = [0.0] * c.maxnbf
    gf[0] = 1.0
    p.TN.beta.fill(0.0)
    for ix in range(c.itb0):
        p.TN.beta[0, ix] = 0.00445 - 0.000022 * ix
    p.TN.beta[0, c.itb0] = 445.0
    p.TN.beta[0, c.itgb0] = 20.0
    p.TN.beta[0, c.itex] = 1200.0
    for s in (p.N2,p.O2,p.O1,p.HE,p.H1,p.AR,p.N1,p.OA,p.NO): s.beta.fill(0.0)
    vals = [
        (p.N2,(100,8,15)), (p.O2,(101,8.5,15.5)), (p.O1,(100,8,15)), (p.HE,(102,9,16)),
        (p.H1,(100,8,16)), (p.AR,(99,8,14)), (p.N1,(101,8,16)), (p.NO,(101,8,16))
    ]
    for s,(zm,hml,hmu) in vals:
        s.beta[0,1]=zm; s.beta[0,2]=hml; s.beta[0,3]=hmu
        s.beta[0,7]=0.03; s.beta[0,8]=106.0; s.beta[0,9]=25.0
    p.O1.beta[0,0]=36.0; p.O1.beta[0,4]=0.15; p.O1.beta[0,5]=92.0; p.O1.beta[0,6]=20.0
    for j in range(c.nsplO1): p.O1.beta[0,j+10]=34.8-0.18*j
    p.H1.beta[0,0]=29.5; p.H1.beta[0,4]=0.12; p.H1.beta[0,5]=95.0; p.H1.beta[0,6]=22.0
    p.N1.beta[0,0]=30.0; p.N1.beta[0,4]=0.11; p.N1.beta[0,5]=96.0; p.N1.beta[0,6]=22.0
    p.OA.beta[0,0]=25.0; p.OA.beta[0,4]=0.08; p.OA.beta[0,5]=180.0; p.OA.beta[0,6]=120.0
    p.NO.beta[0,0]=31.0; p.NO.beta[0,4]=0.10; p.NO.beta[0,5]=96.0; p.NO.beta[0,6]=21.0
    for j in range(c.nsplNO): p.NO.beta[0,j+10]=30.8-0.14*j
    p.N2Rflag=False
    return gf, t.tfnparm(gf)


def _vertical(z,tpro):
    w=[0.0]*4; iz=0
    if z<c.zetaB:
        s,iz=u.bspline(z,c.nodesTN,c.nd,4,p.etaTN); w=[s[-3,4],s[-2,4],s[-1,4],s[0,4]]
    return (t.tfnx(z,iz,w,tpro), tpro.lndtotF-0.0125*(z-c.zetaF),
            tpro.VzetaF+(z-c.zetaF)/max(tpro.tex,300.0),
            tpro.WzetaA+0.001*(z-c.zetaA), 0.80+0.001*min(max(z,0.0),150.0))


class DensityTests(unittest.TestCase):
    def test_dnparm_shape_contract(self):
        x=d.DnParm(); self.assertEqual(len(x.cf),c.nsplO1+2); self.assertEqual(len(x.Mi),5); self.assertEqual(len(x.aMi),5)

    def test_pwmp_endpoint_semantics(self):
        zm=(80.,90.,100.,115.,135.); m=(28.,25.,20.,16.,12.); q=(-.3,-.5,-.2666666666666667,-.2)
        self.assertEqual(d.pwmp(80.,zm,m,q),28.); self.assertEqual(d.pwmp(90.,zm,m,q),25.); self.assertEqual(d.pwmp(135.,zm,m,q),12.)

    def test_invalid_species_rejected(self):
        gf,tp=_seed()
        with self.assertRaises(ValueError): d.dfnparm(1,gf,tp)

    def test_n2rflag_changes_n2_R(self):
        gf,tp=_seed(); p.N2.beta[0,7]=0.25
        p.N2Rflag=False; a=d.dfnparm(2,gf,tp).R
        p.N2Rflag=True; b=d.dfnparm(2,gf,tp).R
        self.assertEqual(a,0.0); self.assertNotEqual(a,b)

    def test_atomic_oxygen_minimum_is_50km(self):
        gf,tp=_seed(); x=d.dfnparm(4,gf,tp); self.assertEqual(x.zmin,50.0)

    def test_no_default_real_cutoff(self):
        gf,tp=_seed(); x=d.dfnparm(10,gf,tp); self.assertEqual(x.zmin,72.5)

    def test_no_missing_branch(self):
        gf,tp=_seed(); p.NO.beta[0,0]=0.0; x=d.dfnparm(10,gf,tp)
        self.assertEqual(x.lndref,0.0); self.assertEqual(d.dfnx(100.0,*_vertical(100.0,tp),tp,x),c.dmissing)

    def test_anomalous_oxygen_cutoff(self):
        gf,tp=_seed(); x=d.dfnparm(9,gf,tp)
        self.assertEqual(d.dfnx(119.999999,*_vertical(119.999999,tp),tp,x),c.dmissing)
        self.assertTrue(math.isfinite(d.dfnx(120.0,*_vertical(120.0,tp),tp,x)))

    def test_effective_mass_nodes_are_ordered(self):
        gf,tp=_seed(); x=d.dfnparm(5,gf,tp); self.assertEqual(x.zetaMi,sorted(x.zetaMi)); self.assertTrue(all(math.isfinite(v) for v in x.Mi+x.aMi+x.WMi+x.XMi))

    def test_oxygen_spline_constraints_present(self):
        gf,tp=_seed(); x=d.dfnparm(4,gf,tp); self.assertNotEqual(x.cf[8],0.0); self.assertNotEqual(x.cf[9],0.0)

    def test_no_spline_constraints_present(self):
        gf,tp=_seed(); x=d.dfnparm(10,gf,tp); self.assertNotEqual(x.cf[8],0.0); self.assertNotEqual(x.cf[9],0.0)

    def test_all_species_representative_density_finite(self):
        gf,tp=_seed()
        for ispec in range(2,11):
            x=d.dfnparm(ispec,gf,tp); z=max(150.0,x.zmin)
            val=d.dfnx(z,*_vertical(z,tp),tp,x)
            self.assertTrue(math.isfinite(val)); self.assertGreaterEqual(val,0.0)

    def test_density_module_has_no_pymsis_dependency(self):
        import inspect
        self.assertNotIn('pymsis', inspect.getsource(d).lower())


if __name__ == '__main__': unittest.main()
