"""
Native NRLMSIS 2.1 vertical temperature formulation.

Authoritative counterpart
-------------------------
NRL NRLMSIS 2.1 ``msis_tfn.F90``.

Derivative translation notice
-----------------------------
This file is a Python translation for the Elara X NRLMSIS native component.
The authoritative source's derived-type fields, logical array indexing,
statement ordering, matrix orientation, spline/Bates boundary semantics, and
dependencies on the frozen Stage-1--4 native layers are preserved.

Use and modification are governed by ``LICENSE_NRLMSIS21.txt`` in the
repository root. See repository provenance and translation-governance
materials for the controlled translation and verification process.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Mapping, MutableSequence, Sequence

from .constants import (
    Mbarg0divkB,
    S4zetaA,
    S4zetaF,
    S5zeta0,
    S5zetaA,
    S5zetaB,
    S5zetaF,
    S6zetaA,
    S6zetaB,
    c2tn,
    cmag,
    cut,
    itb0,
    itex,
    itgb0,
    izax,
    izfx,
    kB,
    lnP0,
    maxnbf,
    mbf,
    nmag,
    nl,
    nut,
    wbeta,
    wgamma,
    wghtAxdz,
    zetaA,
    zetaB,
)
from . import parameters as _parameters
from .horizontal import geomag, sfluxmod, utdep
from .utilities import dilog


@dataclass
class TnParm:
    """Python representation of authoritative Fortran ``type(tnparm)``."""

    cf: list[float] = field(default_factory=lambda: [0.0] * (nl + 1))
    tzetaF: float = 0.0
    tzetaA: float = 0.0
    dlntdzA: float = 0.0
    lndtotF: float = 0.0
    tex: float = 0.0
    tgb0: float = 0.0
    tb0: float = 0.0
    sigma: float = 0.0
    sigmasq: float = 0.0
    b: float = 0.0
    beta: list[float] = field(default_factory=lambda: [0.0] * (nl + 1))
    gamma: list[float] = field(default_factory=lambda: [0.0] * (nl + 1))
    cVs: float = 0.0
    cVb: float = 0.0
    cWs: float = 0.0
    cWb: float = 0.0
    VzetaF: float = 0.0
    VzetaA: float = 0.0
    WzetaA: float = 0.0
    Vzeta0: float = 0.0


# Source-name compatibility for later controlled translations.
tnparm = TnParm


def _require_len(values: Sequence[object], expected: int, name: str) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} must contain exactly {expected} values")


def _dot_beta_gf(col: int, gf: Sequence[float]) -> float:
    beta = _parameters.TN.beta
    if beta is None:
        raise RuntimeError("NRLMSIS parameter space is not initialized")
    total = 0.0
    # Preserve authoritative increasing-index evaluation order.
    for j in range(0, mbf + 1):
        total = total + beta[j, col] * float(gf[j])
    return total


def _dot(values_a: Sequence[float], values_b: Sequence[float]) -> float:
    if len(values_a) != len(values_b):
        raise ValueError("dot-product length mismatch")
    total = 0.0
    for a, b in zip(values_a, values_b):
        total = total + float(a) * float(b)
    return total


def _geomag_plg_from_gf(gf: Sequence[float]) -> list[list[float]]:
    """Map a 14-element Fortran sequence-associated section to (0:6,0:1).

    Fortran storage is column-major, so the first seven elements fill the
    second dummy dimension's first column, followed by the next seven.
    """
    start = cmag + 13
    return [
        [float(gf[start + i]), float(gf[start + 7 + i])]
        for i in range(7)
    ]


def _reset_output(tpro: TnParm) -> None:
    # INTENT(OUT): prior field values are not semantically retained.
    tpro.cf = [0.0] * (nl + 1)
    tpro.beta = [0.0] * (nl + 1)
    tpro.gamma = [0.0] * (nl + 1)
    tpro.tzetaF = 0.0
    tpro.tzetaA = 0.0
    tpro.dlntdzA = 0.0
    tpro.lndtotF = 0.0
    tpro.tex = 0.0
    tpro.tgb0 = 0.0
    tpro.tb0 = 0.0
    tpro.sigma = 0.0
    tpro.sigmasq = 0.0
    tpro.b = 0.0
    tpro.cVs = 0.0
    tpro.cVb = 0.0
    tpro.cWs = 0.0
    tpro.cWb = 0.0
    tpro.VzetaF = 0.0
    tpro.VzetaA = 0.0
    tpro.WzetaA = 0.0
    tpro.Vzeta0 = 0.0


def tfnparm(gf: Sequence[float], tpro: TnParm | None = None) -> TnParm:
    """Compute authoritative vertical temperature/profile parameters."""

    _require_len(gf, maxnbf, "gf")
    if _parameters.TN.beta is None:
        raise RuntimeError("NRLMSIS parameter space is not initialized")

    if tpro is None:
        tpro = TnParm()
    elif not isinstance(tpro, TnParm):
        raise TypeError("tpro must be a TnParm instance")
    _reset_output(tpro)

    # Unconstrained spline coefficients.
    for ix in range(0, itb0):
        tpro.cf[ix] = _dot_beta_gf(ix, gf)

    for ix in range(0, itb0):
        if _parameters.smod[ix]:
            tpro.cf[ix] = tpro.cf[ix] + sfluxmod(
                ix,
                gf,
                _parameters.TN,
                1.0 / _parameters.TN.beta[0, ix],
            )

    plg_input = _geomag_plg_from_gf(gf)

    # Exospheric temperature.
    tpro.tex = _dot_beta_gf(itex, gf)
    tpro.tex = tpro.tex + sfluxmod(
        itex,
        gf,
        _parameters.TN,
        1.0 / _parameters.TN.beta[0, itex],
    )
    tpro.tex = tpro.tex + geomag(
        [_parameters.TN.beta[j, itex] for j in range(cmag, cmag + nmag)],
        [float(gf[j]) for j in range(cmag, cmag + 13)],
        plg_input,
    )
    tpro.tex = tpro.tex + utdep(
        [_parameters.TN.beta[j, itex] for j in range(cut, cut + nut)],
        [float(gf[j]) for j in range(cut, cut + 9)],
    )

    # Temperature gradient at zetaB.
    tpro.tgb0 = _dot_beta_gf(itgb0, gf)
    if _parameters.smod[itgb0]:
        tpro.tgb0 = tpro.tgb0 + sfluxmod(
            itgb0,
            gf,
            _parameters.TN,
            1.0 / _parameters.TN.beta[0, itgb0],
        )
    tpro.tgb0 = tpro.tgb0 + geomag(
        [_parameters.TN.beta[j, itgb0] for j in range(cmag, cmag + nmag)],
        [float(gf[j]) for j in range(cmag, cmag + 13)],
        plg_input,
    )

    # Temperature at zetaB.
    tpro.tb0 = _dot_beta_gf(itb0, gf)
    if _parameters.smod[itb0]:
        tpro.tb0 = tpro.tb0 + sfluxmod(
            itb0,
            gf,
            _parameters.TN,
            1.0 / _parameters.TN.beta[0, itb0],
        )
    tpro.tb0 = tpro.tb0 + geomag(
        [_parameters.TN.beta[j, itb0] for j in range(cmag, cmag + nmag)],
        [float(gf[j]) for j in range(cmag, cmag + 13)],
        plg_input,
    )

    # Shape factor.
    tpro.sigma = tpro.tgb0 / (tpro.tex - tpro.tb0)

    # Constrain the top three spline coefficients for C2 continuity.
    bc = [0.0, 0.0, 0.0]
    bc[0] = 1.0 / tpro.tb0
    bc[1] = -tpro.tgb0 / (tpro.tb0 * tpro.tb0)
    bc[2] = -bc[1] * (tpro.sigma + 2.0 * tpro.tgb0 / tpro.tb0)

    # Fortran MATMUL(bc,c2tn): result(j)=sum_i bc(i)*c2tn(i,j).
    for j in range(3):
        total = 0.0
        for i in range(3):
            total = total + bc[i] * c2tn[i][j]
        tpro.cf[itb0 + j] = total

    # Reference temperature at zetaF.
    tpro.tzetaF = 1.0 / _dot(tpro.cf[izfx : izfx + 3], S4zetaF)

    # Reference temperature and gradient at zetaA.
    tpro.tzetaA = 1.0 / _dot(tpro.cf[izax : izax + 3], S4zetaA)
    tpro.dlntdzA = -_dot(tpro.cf[izax : izax + 3], wghtAxdz) * tpro.tzetaA

    # First and second 1/T integral spline coefficients.
    tpro.beta[0] = tpro.cf[0] * wbeta[0]
    for ix in range(1, nl + 1):
        tpro.beta[ix] = tpro.beta[ix - 1] + tpro.cf[ix] * wbeta[ix]

    tpro.gamma[0] = tpro.beta[0] * wgamma[0]
    for ix in range(1, nl + 1):
        tpro.gamma[ix] = tpro.gamma[ix - 1] + tpro.beta[ix] * wgamma[ix]

    # Integration terms and constants.
    tpro.b = 1.0 - tpro.tb0 / tpro.tex
    tpro.sigmasq = tpro.sigma * tpro.sigma
    tpro.cVs = -_dot(tpro.beta[itb0 - 1 : itb0 + 3], S5zetaB)
    tpro.cWs = -_dot(tpro.gamma[itb0 - 2 : itb0 + 3], S6zetaB)
    tpro.cVb = -math.log(1.0 - tpro.b) / (tpro.sigma * tpro.tex)
    tpro.cWb = -dilog(tpro.b) / (tpro.sigmasq * tpro.tex)
    tpro.VzetaF = _dot(tpro.beta[izfx - 1 : izfx + 3], S5zetaF) + tpro.cVs
    tpro.VzetaA = _dot(tpro.beta[izax - 1 : izax + 3], S5zetaA) + tpro.cVs
    tpro.WzetaA = (
        _dot(tpro.gamma[izax - 2 : izax + 3], S6zetaA)
        + tpro.cVs * (zetaA - zetaB)
        + tpro.cWs
    )
    tpro.Vzeta0 = _dot(tpro.beta[0:3], S5zeta0) + tpro.cVs

    # Total number density at zetaF.
    tpro.lndtotF = (
        lnP0
        - Mbarg0divkB * (tpro.VzetaF - tpro.Vzeta0)
        - math.log(kB * tpro.tzetaF)
    )

    return tpro


def _weight_at(wght: Sequence[float] | Mapping[int, float], logical_index: int) -> float:
    if isinstance(wght, Mapping):
        try:
            return float(wght[logical_index])
        except KeyError as exc:
            raise ValueError("wght mapping must contain logical indices -3..0") from exc
    _require_len(wght, 4, "wght")
    if not -3 <= logical_index <= 0:
        raise ValueError("logical wght index must lie in -3..0")
    return float(wght[logical_index + 3])


def tfnx(
    z: float,
    iz: int,
    wght: Sequence[float] | Mapping[int, float],
    tpro: TnParm,
) -> float:
    """Compute temperature at geopotential height *z* in kilometres."""

    z = float(z)
    iz = int(iz)
    if not isinstance(tpro, TnParm):
        raise TypeError("tpro must be a TnParm instance")

    if z < zetaB:
        if not 0 <= iz <= nl:
            raise ValueError("iz must lie in 0..nl in the spline region")
        i = max(iz - 3, 0)
        if iz < 3:
            j = -iz
        else:
            j = -3

        total = 0.0
        logical_weight = j
        for cf_index in range(i, iz + 1):
            total = total + tpro.cf[cf_index] * _weight_at(wght, logical_weight)
            logical_weight += 1
        return 1.0 / total

    # Exactly zetaB belongs to the Bates region, as in authoritative ".lt.".
    return tpro.tex - (tpro.tex - tpro.tb0) * math.exp(
        -tpro.sigma * (z - zetaB)
    )


__all__ = [
    "TnParm",
    "tnparm",
    "tfnparm",
    "tfnx",
]
