"""
Native NRLMSIS 2.1 vertical species-density formulation.

Authoritative counterpart
-------------------------
NRL NRLMSIS 2.1 ``msis_dfn.F90``.

Derivative translation notice
-----------------------------
This file is a Python translation for the Elara X NRLMSIS native component.
The authoritative species dispatch, derived-type fields, logical indexing,
piecewise effective-mass profile, spline/hydrostatic boundary semantics,
statement ordering, and dependencies on the frozen Stage-1--5 native layers
are preserved.

Use and modification are governed by ``LICENSE_NRLMSIS21.txt`` in the
repository root. See repository provenance and translation-governance
materials for the controlled translation and verification process.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import struct
from typing import Sequence

from .constants import (
    HOA,
    Mbar,
    c1NO,
    c1NOadj,
    c1o1,
    c1o1adj,
    cmag,
    cut,
    dmissing,
    g0divkB,
    lnvmr,
    maxnbf,
    mbf,
    nd,
    ndNO,
    ndO1,
    nmag,
    nodesNO,
    nodesO1,
    nodesTN,
    nsplNO,
    nsplO1,
    nut,
    specmass,
    tanh1,
    zetaA,
    zetaB,
    zetaF,
    zetarefNO,
    zetarefO1,
    zetarefOA,
)
from . import parameters as _parameters
from .horizontal import geomag, sfluxmod, utdep
from .temperature import TnParm
from .utilities import bspline, dilog


def _f32(value: float) -> float:
    """Round through authoritative default REAL before binary64 promotion."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


@dataclass
class DnParm:
    """Python representation of authoritative Fortran ``type(dnparm)``."""

    lnPhiF: float = 0.0
    lndref: float = 0.0
    zetaM: float = 0.0
    HML: float = 0.0
    HMU: float = 0.0
    C: float = 0.0
    zetaC: float = 0.0
    HC: float = 0.0
    R: float = 0.0
    zetaR: float = 0.0
    HR: float = 0.0
    cf: list[float] = field(default_factory=lambda: [0.0] * (nsplO1 + 2))
    zref: float = 0.0
    Mi: list[float] = field(default_factory=lambda: [0.0] * 5)
    zetaMi: list[float] = field(default_factory=lambda: [0.0] * 5)
    aMi: list[float] = field(default_factory=lambda: [0.0] * 5)
    WMi: list[float] = field(default_factory=lambda: [0.0] * 5)
    XMi: list[float] = field(default_factory=lambda: [0.0] * 5)
    Izref: float = 0.0
    Tref: float = 0.0
    zmin: float = 0.0
    zhyd: float = 0.0
    ispec: int = 0


# Source-name compatibility for later controlled translations.
dnparm = DnParm


def _require_len(values: Sequence[object], expected: int, name: str) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} must contain exactly {expected} values")


def _reset_output(dpro: DnParm) -> None:
    # INTENT(OUT): prior values are not semantically retained.
    dpro.lnPhiF = 0.0
    dpro.lndref = 0.0
    dpro.zetaM = 0.0
    dpro.HML = 0.0
    dpro.HMU = 0.0
    dpro.C = 0.0
    dpro.zetaC = 0.0
    dpro.HC = 0.0
    dpro.R = 0.0
    dpro.zetaR = 0.0
    dpro.HR = 0.0
    dpro.cf = [0.0] * (nsplO1 + 2)
    dpro.zref = 0.0
    dpro.Mi = [0.0] * 5
    dpro.zetaMi = [0.0] * 5
    dpro.aMi = [0.0] * 5
    dpro.WMi = [0.0] * 5
    dpro.XMi = [0.0] * 5
    dpro.Izref = 0.0
    dpro.Tref = 0.0
    dpro.zmin = 0.0
    dpro.zhyd = 0.0
    dpro.ispec = 0


def _dot_subset(subset, col: int, gf: Sequence[float]) -> float:
    if subset.beta is None:
        raise RuntimeError("NRLMSIS parameter space is not initialized")
    total = 0.0
    for j in range(0, mbf + 1):
        total = total + subset.beta[j, col] * float(gf[j])
    return total


def _geomag_plg_from_gf(gf: Sequence[float]) -> list[list[float]]:
    start = cmag + 13
    return [
        [float(gf[start + i]), float(gf[start + 7 + i])]
        for i in range(7)
    ]


def _geomag_subset(subset, col: int, gf: Sequence[float]) -> float:
    if subset.beta is None:
        raise RuntimeError("NRLMSIS parameter space is not initialized")
    return geomag(
        [subset.beta[j, col] for j in range(cmag, cmag + nmag)],
        [float(gf[j]) for j in range(cmag, cmag + 13)],
        _geomag_plg_from_gf(gf),
    )


def _ut_subset(subset, col: int, gf: Sequence[float]) -> float:
    if subset.beta is None:
        raise RuntimeError("NRLMSIS parameter space is not initialized")
    return utdep(
        [subset.beta[j, col] for j in range(cut, cut + nut)],
        [float(gf[j]) for j in range(cut, cut + 9)],
    )


def _matmul_row2(vec: Sequence[float], matrix: Sequence[Sequence[float]]) -> list[float]:
    if len(vec) != 2 or len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError("2x2 row-vector MATMUL shape mismatch")
    result = [0.0, 0.0]
    for j in range(2):
        total = 0.0
        for i in range(2):
            total = total + float(vec[i]) * float(matrix[i][j])
        result[j] = total
    return result


def pwmp(z: float, zm: Sequence[float], m: Sequence[float], dmdz: Sequence[float]) -> float:
    """Piecewise effective-mass interpolation with Fortran equality semantics."""
    _require_len(zm, 5, "zm")
    _require_len(m, 5, "m")
    if len(dmdz) < 4:
        raise ValueError("dmdz must contain at least four values")
    z = float(z)
    if z >= float(zm[4]):
        return float(m[4])
    if z <= float(zm[0]):
        return float(m[0])
    for inode in range(0, 4):
        if z < float(zm[inode + 1]):
            return float(m[inode]) + float(dmdz[inode]) * (z - float(zm[inode]))
    raise RuntimeError("authoritative PWMP interval resolution failed")


def dfnparm(
    ispec: int,
    gf: Sequence[float],
    tpro: TnParm,
    dpro: DnParm | None = None,
) -> DnParm:
    """Compute authoritative vertical species-density profile parameters."""
    _require_len(gf, maxnbf, "gf")
    if not isinstance(tpro, TnParm):
        raise TypeError("tpro must be a TnParm instance")
    ispec = int(ispec)
    if ispec < 2 or ispec > 10:
        raise ValueError("species index must lie in 2..10")
    if dpro is None:
        dpro = DnParm()
    elif not isinstance(dpro, DnParm):
        raise TypeError("dpro must be a DnParm instance")
    _reset_output(dpro)
    dpro.ispec = ispec

    p = _parameters

    if ispec == 2:  # N2
        dpro.lnPhiF = lnvmr[ispec - 1]
        dpro.lndref = tpro.lndtotF + dpro.lnPhiF
        dpro.zref = zetaF
        dpro.zmin = -1.0
        dpro.zhyd = zetaF
        dpro.zetaM = _dot_subset(p.N2, 1, gf)
        dpro.HML = p.N2.beta[0, 2]
        dpro.HMU = p.N2.beta[0, 3]
        dpro.R = 0.0
        if p.N2Rflag:
            dpro.R = _dot_subset(p.N2, 7, gf)
        dpro.zetaR = p.N2.beta[0, 8]
        dpro.HR = p.N2.beta[0, 9]

    elif ispec == 3:  # O2
        dpro.lnPhiF = lnvmr[ispec - 1]
        dpro.lndref = tpro.lndtotF + dpro.lnPhiF
        dpro.zref = zetaF
        dpro.zmin = -1.0
        dpro.zhyd = zetaF
        dpro.zetaM = p.O2.beta[0, 1]
        dpro.HML = p.O2.beta[0, 2]
        dpro.HMU = p.O2.beta[0, 3]
        dpro.R = _dot_subset(p.O2, 7, gf)
        dpro.R = dpro.R + _geomag_subset(p.O2, 7, gf)
        dpro.zetaR = p.O2.beta[0, 8]
        dpro.HR = p.O2.beta[0, 9]

    elif ispec == 4:  # O
        dpro.lnPhiF = 0.0
        dpro.lndref = _dot_subset(p.O1, 0, gf)
        dpro.zref = zetarefO1
        dpro.zmin = nodesO1[3]
        dpro.zhyd = zetarefO1
        dpro.zetaM = p.O1.beta[0, 1]
        dpro.HML = p.O1.beta[0, 2]
        dpro.HMU = p.O1.beta[0, 3]
        dpro.C = _dot_subset(p.O1, 4, gf)
        dpro.zetaC = p.O1.beta[0, 5]
        dpro.HC = p.O1.beta[0, 6]
        dpro.R = _dot_subset(p.O1, 7, gf)
        dpro.R = dpro.R + sfluxmod(7, gf, p.O1, 0.0)
        dpro.R = dpro.R + _geomag_subset(p.O1, 7, gf)
        dpro.R = dpro.R + _ut_subset(p.O1, 7, gf)
        dpro.zetaR = p.O1.beta[0, 8]
        dpro.HR = p.O1.beta[0, 9]
        for izf in range(0, nsplO1):
            dpro.cf[izf] = _dot_subset(p.O1, izf + 10, gf)

    elif ispec == 5:  # He
        dpro.lnPhiF = lnvmr[ispec - 1]
        dpro.lndref = tpro.lndtotF + dpro.lnPhiF
        dpro.zref = zetaF
        dpro.zmin = -1.0
        dpro.zhyd = zetaF
        dpro.zetaM = p.HE.beta[0, 1]
        dpro.HML = p.HE.beta[0, 2]
        dpro.HMU = p.HE.beta[0, 3]
        dpro.R = _dot_subset(p.HE, 7, gf)
        dpro.R = dpro.R + sfluxmod(7, gf, p.HE, 1.0)
        dpro.R = dpro.R + _geomag_subset(p.HE, 7, gf)
        dpro.R = dpro.R + _ut_subset(p.HE, 7, gf)
        dpro.zetaR = p.HE.beta[0, 8]
        dpro.HR = p.HE.beta[0, 9]

    elif ispec == 6:  # H
        dpro.lnPhiF = 0.0
        dpro.lndref = _dot_subset(p.H1, 0, gf)
        dpro.zref = zetaA
        dpro.zmin = 75.0
        dpro.zhyd = zetaF
        dpro.zetaM = p.H1.beta[0, 1]
        dpro.HML = p.H1.beta[0, 2]
        dpro.HMU = p.H1.beta[0, 3]
        dpro.C = _dot_subset(p.H1, 4, gf)
        dpro.zetaC = _dot_subset(p.H1, 5, gf)
        dpro.HC = p.H1.beta[0, 6]
        dpro.R = _dot_subset(p.H1, 7, gf)
        dpro.R = dpro.R + sfluxmod(7, gf, p.H1, 0.0)
        dpro.R = dpro.R + _geomag_subset(p.H1, 7, gf)
        dpro.R = dpro.R + _ut_subset(p.H1, 7, gf)
        dpro.zetaR = p.H1.beta[0, 8]
        dpro.HR = p.H1.beta[0, 9]

    elif ispec == 7:  # Ar
        dpro.lnPhiF = lnvmr[ispec - 1]
        dpro.lndref = tpro.lndtotF + dpro.lnPhiF
        dpro.zref = zetaF
        dpro.zmin = -1.0
        dpro.zhyd = zetaF
        dpro.zetaM = p.AR.beta[0, 1]
        dpro.HML = p.AR.beta[0, 2]
        dpro.HMU = p.AR.beta[0, 3]
        dpro.R = _dot_subset(p.AR, 7, gf)
        dpro.R = dpro.R + _geomag_subset(p.AR, 7, gf)
        dpro.R = dpro.R + _ut_subset(p.AR, 7, gf)
        dpro.zetaR = p.AR.beta[0, 8]
        dpro.HR = p.AR.beta[0, 9]

    elif ispec == 8:  # N
        dpro.lnPhiF = 0.0
        dpro.lndref = _dot_subset(p.N1, 0, gf)
        dpro.lndref = dpro.lndref + sfluxmod(0, gf, p.N1, 0.0)
        dpro.lndref = dpro.lndref + _geomag_subset(p.N1, 0, gf)
        dpro.lndref = dpro.lndref + _ut_subset(p.N1, 0, gf)
        dpro.zref = zetaB
        dpro.zmin = 90.0
        dpro.zhyd = zetaF
        dpro.zetaM = p.N1.beta[0, 1]
        dpro.HML = p.N1.beta[0, 2]
        dpro.HMU = p.N1.beta[0, 3]
        dpro.C = p.N1.beta[0, 4]
        dpro.zetaC = p.N1.beta[0, 5]
        dpro.HC = p.N1.beta[0, 6]
        dpro.R = _dot_subset(p.N1, 7, gf)
        dpro.zetaR = p.N1.beta[0, 8]
        dpro.HR = p.N1.beta[0, 9]

    elif ispec == 9:  # anomalous O
        dpro.lndref = _dot_subset(p.OA, 0, gf)
        dpro.lndref = dpro.lndref + _geomag_subset(p.OA, 0, gf)
        dpro.zref = zetarefOA
        dpro.zmin = 120.0
        dpro.zhyd = 0.0
        dpro.C = p.OA.beta[0, 4]
        dpro.zetaC = p.OA.beta[0, 5]
        dpro.HC = p.OA.beta[0, 6]
        return dpro

    elif ispec == 10:  # NO
        if p.NO.beta[0, 0] == 0.0:
            dpro.lndref = 0.0
            return dpro
        dpro.lnPhiF = 0.0
        dpro.lndref = _dot_subset(p.NO, 0, gf)
        dpro.lndref = dpro.lndref + _geomag_subset(p.NO, 0, gf)
        dpro.zref = zetarefNO
        # Authoritative literal is unsuffixed default REAL before promotion.
        dpro.zmin = _f32(72.5)
        dpro.zhyd = zetarefNO
        dpro.zetaM = _dot_subset(p.NO, 1, gf)
        dpro.HML = _dot_subset(p.NO, 2, gf)
        dpro.HMU = _dot_subset(p.NO, 3, gf)
        dpro.C = _dot_subset(p.NO, 4, gf)
        dpro.C = dpro.C + _geomag_subset(p.NO, 4, gf)
        dpro.zetaC = _dot_subset(p.NO, 5, gf)
        dpro.HC = _dot_subset(p.NO, 6, gf)
        dpro.R = _dot_subset(p.NO, 7, gf)
        dpro.zetaR = _dot_subset(p.NO, 8, gf)
        dpro.HR = _dot_subset(p.NO, 9, gf)
        for izf in range(0, nsplNO):
            dpro.cf[izf] = _dot_subset(p.NO, izf + 10, gf)
            dpro.cf[izf] = dpro.cf[izf] + _geomag_subset(p.NO, izf + 10, gf)

    # Piecewise effective-mass profile and integration terms.
    dpro.zetaMi[0] = dpro.zetaM - 2.0 * dpro.HML
    dpro.zetaMi[1] = dpro.zetaM - dpro.HML
    dpro.zetaMi[2] = dpro.zetaM
    dpro.zetaMi[3] = dpro.zetaM + dpro.HMU
    dpro.zetaMi[4] = dpro.zetaM + 2.0 * dpro.HMU
    dpro.Mi[0] = Mbar
    dpro.Mi[4] = specmass[ispec - 1]
    dpro.Mi[2] = (dpro.Mi[0] + dpro.Mi[4]) / 2.0
    delM = tanh1 * (dpro.Mi[4] - dpro.Mi[0]) / 2.0
    dpro.Mi[1] = dpro.Mi[2] - delM
    dpro.Mi[3] = dpro.Mi[2] + delM
    for i in range(0, 4):
        dpro.aMi[i] = (dpro.Mi[i + 1] - dpro.Mi[i]) / (dpro.zetaMi[i + 1] - dpro.zetaMi[i])

    for i in range(0, 5):
        delz = dpro.zetaMi[i] - zetaB
        if dpro.zetaMi[i] < zetaB:
            s, iz = bspline(dpro.zetaMi[i], nodesTN, nd + 2, 6, p.etaTN)
            total = 0.0
            for logical in range(-5, 1):
                total = total + tpro.gamma[iz + logical] * s[logical, 6]
            dpro.WMi[i] = total + tpro.cVs * delz + tpro.cWs
        else:
            dpro.WMi[i] = (
                0.5 * delz * delz
                + dilog(tpro.b * math.exp(-tpro.sigma * delz)) / tpro.sigmasq
            ) / tpro.tex + tpro.cVb * delz + tpro.cWb

    dpro.XMi[0] = -dpro.aMi[0] * dpro.WMi[0]
    for i in range(1, 4):
        dpro.XMi[i] = dpro.XMi[i - 1] - dpro.WMi[i] * (dpro.aMi[i] - dpro.aMi[i - 1])
    dpro.XMi[4] = dpro.XMi[3] + dpro.WMi[4] * dpro.aMi[3]

    # Hydrostatic integral and reference temperature.
    if dpro.zref == zetaF:
        Mzref = Mbar
        dpro.Tref = tpro.tzetaF
        dpro.Izref = Mbar * tpro.VzetaF
    elif dpro.zref == zetaB:
        Mzref = pwmp(dpro.zref, dpro.zetaMi, dpro.Mi, dpro.aMi)
        dpro.Tref = tpro.tb0
        dpro.Izref = 0.0
        if zetaB > dpro.zetaMi[0] and zetaB < dpro.zetaMi[4]:
            i = 0
            for i1 in range(1, 4):
                if zetaB < dpro.zetaMi[i1]:
                    break
                i = i1
            dpro.Izref = dpro.Izref - dpro.XMi[i]
        else:
            dpro.Izref = dpro.Izref - dpro.XMi[4]
    elif dpro.zref == zetaA:
        Mzref = pwmp(dpro.zref, dpro.zetaMi, dpro.Mi, dpro.aMi)
        dpro.Tref = tpro.tzetaA
        dpro.Izref = Mzref * tpro.VzetaA
        if zetaA > dpro.zetaMi[0] and zetaA < dpro.zetaMi[4]:
            i = 0
            for i1 in range(1, 4):
                if zetaA < dpro.zetaMi[i1]:
                    break
                i = i1
            dpro.Izref = dpro.Izref - (dpro.aMi[i] * tpro.WzetaA + dpro.XMi[i])
        else:
            dpro.Izref = dpro.Izref - dpro.XMi[4]
    else:
        raise RuntimeError("authoritative hydrostatic integral unavailable at reference height")

    # C1 constraint for O at 85 km.
    if ispec == 4:
        Cterm = dpro.C * math.exp(-(dpro.zref - dpro.zetaC) / dpro.HC)
        Rterm0 = math.tanh((dpro.zref - dpro.zetaR) / (p.HRfactO1ref * dpro.HR))
        Rterm = dpro.R * (1.0 + Rterm0)
        bc0 = dpro.lndref - Cterm + Rterm - dpro.cf[7] * c1o1adj[0]
        bc1 = (
            -Mzref * g0divkB / tpro.tzetaA
            - tpro.dlntdzA
            + Cterm / dpro.HC
            + Rterm * (1.0 - Rterm0) / dpro.HR * p.dHRfactO1ref
            - dpro.cf[7] * c1o1adj[1]
        )
        dpro.cf[8:10] = _matmul_row2((bc0, bc1), c1o1)

    # C1 constraint for NO at 122.5 km.
    if ispec == 10:
        Cterm = dpro.C * math.exp(-(dpro.zref - dpro.zetaC) / dpro.HC)
        Rterm0 = math.tanh((dpro.zref - dpro.zetaR) / (p.HRfactNOref * dpro.HR))
        Rterm = dpro.R * (1.0 + Rterm0)
        bc0 = dpro.lndref - Cterm + Rterm - dpro.cf[7] * c1NOadj[0]
        bc1 = (
            -Mzref * g0divkB / tpro.tb0
            - tpro.tgb0 / tpro.tb0
            + Cterm / dpro.HC
            + Rterm * (1.0 - Rterm0) / dpro.HR * p.dHRfactNOref
            - dpro.cf[7] * c1NOadj[1]
        )
        dpro.cf[8:10] = _matmul_row2((bc0, bc1), c1NO)

    return dpro


def dfnx(
    z: float,
    tnz: float,
    lndtotz: float,
    Vz: float,
    Wz: float,
    HRfact: float,
    tpro: TnParm,
    dpro: DnParm,
) -> float:
    """Compute species number density at geopotential height *z*."""
    if not isinstance(tpro, TnParm):
        raise TypeError("tpro must be a TnParm instance")
    if not isinstance(dpro, DnParm):
        raise TypeError("dpro must be a DnParm instance")

    z = float(z)
    tnz = float(tnz)
    lndtotz = float(lndtotz)
    Vz = float(Vz)
    Wz = float(Wz)
    HRfact = float(HRfact)

    if z < dpro.zmin:
        return dmissing

    if dpro.ispec == 9:
        value = dpro.lndref - (z - dpro.zref) / HOA - dpro.C * math.exp(-(z - dpro.zetaC) / dpro.HC)
        return math.exp(value)

    if dpro.ispec == 10 and dpro.lndref == 0.0:
        return dmissing

    if dpro.ispec in (2, 3, 5, 7):
        ccor = dpro.R * (1.0 + math.tanh((z - dpro.zetaR) / (HRfact * dpro.HR)))
    elif dpro.ispec in (4, 6, 8, 10):
        ccor = (
            -dpro.C * math.exp(-(z - dpro.zetaC) / dpro.HC)
            + dpro.R * (1.0 + math.tanh((z - dpro.zetaR) / (HRfact * dpro.HR)))
        )
    else:
        raise ValueError("unsupported species index in dfnx")

    if z < dpro.zhyd:
        if dpro.ispec in (2, 3, 5, 7):
            return math.exp(lndtotz + dpro.lnPhiF + ccor)
        if dpro.ispec == 4:
            s, iz = bspline(z, nodesO1, ndO1, 4, _parameters.etaO1)
            total = 0.0
            logical = -3
            for cf_index in range(iz - 3, iz + 1):
                total = total + dpro.cf[cf_index] * s[logical, 4]
                logical += 1
            return math.exp(total)
        if dpro.ispec == 10:
            s, iz = bspline(z, nodesNO, ndNO, 4, _parameters.etaNO)
            total = 0.0
            logical = -3
            for cf_index in range(iz - 3, iz + 1):
                total = total + dpro.cf[cf_index] * s[logical, 4]
                logical += 1
            return math.exp(total)

    Mz = pwmp(z, dpro.zetaMi, dpro.Mi, dpro.aMi)
    Ihyd = Mz * Vz - dpro.Izref
    if z > dpro.zetaMi[0] and z < dpro.zetaMi[4]:
        i = 0
        for i1 in range(1, 4):
            if z < dpro.zetaMi[i1]:
                break
            i = i1
        Ihyd = Ihyd - (dpro.aMi[i] * Wz + dpro.XMi[i])
    elif z >= dpro.zetaMi[4]:
        Ihyd = Ihyd - dpro.XMi[4]

    value = dpro.lndref - Ihyd * g0divkB + ccor
    return math.exp(value) * dpro.Tref / tnz


__all__ = [
    "DnParm",
    "dnparm",
    "dfnparm",
    "dfnx",
    "pwmp",
]
