"""
Native NRLMSIS 2.1 main model-calculation orchestration.

Authoritative counterpart
-------------------------
NRL NRLMSIS 2.1 ``msis_calc.F90``.

Derivative translation notice
-----------------------------
This file is a Python translation for the Elara X NRLMSIS native component.
The authoritative input/output ordering, altitude conversion, SAVE/cache
semantics, spline/hydrostatic branch boundaries, species gating, mass-density
calculation and optional exospheric-temperature semantics are preserved.

Parameter-resource acquisition/distribution remains a separate controlled
resource-layer concern.  As in the authoritative source, ``msiscalc`` invokes
``parameters.msisinit()`` only when the model has not already been initialized.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import struct
from typing import Sequence

from .constants import (
    Hgamma,
    Mbarg0divkB,
    dmissing,
    kB,
    lnP0,
    maxnbf,
    nd,
    nodesTN,
    nspec,
    zetaB,
    zetaF,
    zetagamma,
)
from . import parameters as _parameters
from .density import DnParm, dfnparm, dfnx
from .horizontal import globe
from .temperature import TnParm, tfnparm, tfnx
from .utilities import alt2gph, bspline, dilog


def _f32(value: float) -> float:
    """Round a source default-REAL literal before binary64 promotion."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def _dot(values_a: Sequence[float], values_b: Sequence[float]) -> float:
    """Fortran-like left-to-right scalar dot product for contiguous slices."""
    if len(values_a) != len(values_b):
        raise ValueError("dot-product length mismatch")
    total = 0.0
    for a, b in zip(values_a, values_b):
        total = total + float(a) * float(b)
    return total


def _require_ap(ap: Sequence[float]) -> tuple[float, ...]:
    if len(ap) != 7:
        raise ValueError("ap must contain exactly 7 values")
    return tuple(float(x) for x in ap)


@dataclass
class _ModelCache:
    """Python representation of MSISCALC local SAVE state."""

    lastday: float = field(default_factory=lambda: _f32(-9999.0))
    lastutsec: float = field(default_factory=lambda: _f32(-9999.0))
    lastlat: float = field(default_factory=lambda: _f32(-9999.0))
    lastlon: float = field(default_factory=lambda: _f32(-9999.0))
    lastz: float = field(default_factory=lambda: _f32(-9999.0))
    lastsflux: float = field(default_factory=lambda: _f32(-9999.0))
    lastsfluxavg: float = field(default_factory=lambda: _f32(-9999.0))
    lastap: list[float] = field(default_factory=lambda: [_f32(-9999.0)] * 7)
    gf: list[float] = field(default_factory=lambda: [0.0] * maxnbf)
    Sz: dict[tuple[int, int], float] = field(
        default_factory=lambda: {(l, k): 0.0 for k in range(2, 7) for l in range(-5, 1)}
    )
    iz: int = 0
    tpro: TnParm = field(default_factory=TnParm)
    dpro: list[DnParm | None] = field(
        default_factory=lambda: [None] + [DnParm() for _ in range(nspec - 1)]
    )


_CACHE = _ModelCache()


def _reset_cache_for_testing() -> None:
    """Reset MSISCALC SAVE state for controlled tests/oracle generation.

    This is intentionally private and is not part of the scientific public API.
    """
    global _CACHE
    _CACHE = _ModelCache()


def _profile_changed(
    day: float,
    utsec: float,
    lat: float,
    lon: float,
    sfluxavg: float,
    sflux: float,
    ap: Sequence[float],
) -> bool:
    c = _CACHE
    return (
        day != c.lastday
        or utsec != c.lastutsec
        or lat != c.lastlat
        or lon != c.lastlon
        or sflux != c.lastsflux
        or sfluxavg != c.lastsfluxavg
        or any(a != b for a, b in zip(ap, c.lastap))
    )


def msiscalc(
    day: float,
    utsec: float,
    z: float,
    lat: float,
    lon: float,
    sfluxavg: float,
    sflux: float,
    ap: Sequence[float],
    *,
    return_tex: bool = False,
):
    """Evaluate the native NRLMSIS 2.1 main calculation.

    Parameters follow authoritative ``MSISCALC`` ordering.  ``ap`` contains the
    seven authoritative Ap-history values.  The native Python return is
    ``(tn, dn)`` by default and ``(tn, dn, tex)`` when ``return_tex=True``;
    ``dn`` is a 10-element tuple corresponding exactly to Fortran ``DN(1:10)``.
    """

    day = float(day)
    utsec = float(utsec)
    z = float(z)
    lat = float(lat)
    lon = float(lon)
    sfluxavg = float(sfluxavg)
    sflux = float(sflux)
    apv = _require_ap(ap)

    # Authoritative MSISCALC self-initialization contract.  Stage-7 controlled
    # equivalence uses preinitialized synthetic state; resource policy is Stage 9.
    if not _parameters.initflag:
        _parameters.msisinit()

    # In the DBLE authoritative build, dble(z/lat) feeds the REAL(8) ALT2GPH
    # utility.  Python float is binary64, so float conversion preserves this path.
    if _parameters.zaltflag:
        zeta = alt2gph(float(lat), float(z))
    else:
        zeta = z

    c = _CACHE

    # Altitude-only SAVE cache: exactly zetaB belongs to the analytical branch.
    if zeta < zetaB:
        if zeta != c.lastz:
            kmax = 5 if zeta < zetaF else 6
            c.Sz, c.iz = bspline(zeta, nodesTN, nd + 2, kmax, _parameters.etaTN)
            c.lastz = zeta

    # Profile cache uses exact Fortran .NE. semantics for each scalar and ANY(AP).
    if _profile_changed(day, utsec, lat, lon, sfluxavg, sflux, apv):
        globe(day, utsec, lat, lon, sfluxavg, sflux, apv, c.gf)
        tfnparm(c.gf, c.tpro)
        for ispec in range(2, nspec):
            # Fortran specflag(1:10) maps to Python specflag[0:10].
            if _parameters.specflag[ispec - 1]:
                dpro = c.dpro[ispec]
                if dpro is None:
                    dpro = DnParm()
                    c.dpro[ispec] = dpro
                dfnparm(ispec, c.gf, c.tpro, dpro)

        c.lastday = day
        c.lastutsec = utsec
        c.lastlat = lat
        c.lastlon = lon
        c.lastsflux = sflux
        c.lastsfluxavg = sfluxavg
        c.lastap[:] = apv

    tex = c.tpro.tex

    # Temperature at altitude.  The spline weights carry Fortran logical -3:0.
    w4 = {logical: c.Sz[logical, 4] for logical in range(-3, 1)}
    tn = tfnx(zeta, c.iz, w4, c.tpro)

    # Temperature integration terms. lndtotz is only consumed by density paths
    # below zetaF; initialize it deterministically for Python even though the
    # authoritative local is undefined on upper branches where it is unused.
    lndtotz = 0.0
    delz = zeta - zetaB
    if zeta < zetaF:
        i = max(c.iz - 4, 0)
        j = -c.iz if c.iz < 4 else -4
        beta_values = c.tpro.beta[i : c.iz + 1]
        spline_values = [c.Sz[logical, 5] for logical in range(j, 1)]
        Vz = _dot(beta_values, spline_values) + c.tpro.cVs
        Wz = 0.0
        lnPz = lnP0 - Mbarg0divkB * (Vz - c.tpro.Vzeta0)
        lndtotz = lnPz - math.log(kB * tn)
    else:
        if zeta < zetaB:
            Vz = _dot(
                c.tpro.beta[c.iz - 4 : c.iz + 1],
                [c.Sz[logical, 5] for logical in range(-4, 1)],
            ) + c.tpro.cVs
            Wz = (
                _dot(
                    c.tpro.gamma[c.iz - 5 : c.iz + 1],
                    [c.Sz[logical, 6] for logical in range(-5, 1)],
                )
                + c.tpro.cVs * delz
                + c.tpro.cWs
            )
        else:
            Vz = (
                delz + math.log(tn / c.tpro.tex) / c.tpro.sigma
            ) / c.tpro.tex + c.tpro.cVb
            Wz = (
                0.5 * delz * delz
                + dilog(c.tpro.b * math.exp(-c.tpro.sigma * delz)) / c.tpro.sigmasq
            ) / c.tpro.tex + c.tpro.cVb * delz + c.tpro.cWb

    HRfact = 0.5 * (1.0 + math.tanh(Hgamma * (zeta - zetagamma)))

    dn = [0.0] * 10
    for ispec in range(2, nspec):
        out_index = ispec - 1
        if _parameters.specflag[out_index]:
            dpro = c.dpro[ispec]
            if dpro is None:
                raise RuntimeError("density profile cache missing for enabled species")
            dn[out_index] = dfnx(
                zeta,
                tn,
                lndtotz,
                Vz,
                Wz,
                HRfact,
                c.tpro,
                dpro,
            )
        else:
            dn[out_index] = dmissing

    if _parameters.specflag[0]:
        dn[0] = _dot(dn, _parameters.masswgt)
    else:
        dn[0] = dmissing

    result_dn = tuple(dn)
    if return_tex:
        return tn, result_dn, tex
    return tn, result_dn


__all__ = ["msiscalc"]
