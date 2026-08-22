"""
Native NRLMSIS 2.1 horizontal and time-dependent basis functions.

Authoritative counterpart
-------------------------
NRL NRLMSIS 2.1 ``msis_gfn.F90``.

Derivative translation notice
-----------------------------
This file is a Python translation for the Elara X NRLMSIS native component.
The authoritative source's logical indexing, switch masks, persistent cache
state, evaluation order, and default-REAL promotion semantics are preserved
where they are scientifically relevant.

Use and modification are governed by ``LICENSE_NRLMSIS21.txt`` in the
repository root. See repository provenance and translation-governance
materials for the controlled translation and verification process.
"""

from __future__ import annotations

import math
import struct
from typing import MutableSequence, Sequence

from .constants import (
    amaxn,
    amaxs,
    cextra,
    cintann,
    cmag,
    cnonlin,
    cspw,
    csfx,
    csfxmod,
    ctide,
    ctimeind,
    cut,
    deg2rad,
    doy2rad,
    lst2rad,
    maxn,
    maxnbf,
    mbf,
    nmag,
    nsfx,
    nsfxmod,
    nut,
    pi,
    pmaxm,
    pmaxn,
    pmaxs,
    tmaxl,
    tmaxn,
    tmaxs,
)
from . import parameters as _parameters


def _f32(value: float) -> float:
    """Round a source default-REAL literal to IEEE-754 binary32."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


# Fortran module arrays retain their logical lower bounds by reserving index 0
# where the source starts at 1. PLG is naturally zero-based in both dimensions.
plg = [[0.0 for _m in range(maxn + 1)] for _n in range(maxn + 1)]
cdoy = [0.0, 0.0, 0.0]
sdoy = [0.0, 0.0, 0.0]
clst = [0.0, 0.0, 0.0, 0.0]
slst = [0.0, 0.0, 0.0, 0.0]
clon = [0.0, 0.0, 0.0]
slon = [0.0, 0.0, 0.0]

# These source literals are unsuffixed default REAL values before promotion.
sfluxavgref = _f32(150.0)
sfluxavg_quad_cutoff = _f32(150.0)
lastlat = _f32(-999.9)
lastdoy = _f32(-999.9)
lastlst = _f32(-999.9)
lastlon = _f32(-999.9)

# SOLZEN's unsuffixed parameter-array literals are likewise rounded to default
# REAL before assignment to real(kind=rp) in a DBLE build.
_SOLZEN_P = (
    _f32(0.017203534),
    _f32(0.034407068),
    _f32(0.051610602),
    _f32(0.068814136),
    _f32(0.103221204),
)


def _require_len(values: Sequence[object], expected: int, name: str) -> None:
    if len(values) != expected:
        raise ValueError(f"{name} must contain exactly {expected} values")


def _fortran_mod(a: float, p: float) -> float:
    """Fortran MOD semantics for real arguments (truncation toward zero)."""
    return math.fmod(a, p)


def globe(
    doy: float,
    utsec: float,
    lat: float,
    lon: float,
    sfluxavg: float,
    sflux: float,
    ap: Sequence[float],
    bf: MutableSequence[float] | None = None,
) -> MutableSequence[float]:
    """Calculate the 512 horizontal/time-dependent basis-function terms.

    If *bf* is supplied, it is filled in-place to mirror the Fortran
    ``INTENT(OUT)`` argument. Otherwise a new 512-element list is returned.
    """
    global lastlat, lastdoy, lastlst, lastlon

    _require_len(ap, 7, "ap")
    if bf is None:
        bf = [0.0] * maxnbf
    else:
        _require_len(bf, maxnbf, "bf")
        for j in range(maxnbf):
            bf[j] = 0.0

    doy = float(doy)
    utsec = float(utsec)
    lat = float(lat)
    lon = float(lon)
    sfluxavg = float(sfluxavg)
    sflux = float(sflux)

    if lat != lastlat:
        clat = math.sin(lat * deg2rad)
        slat = math.cos(lat * deg2rad)
        clat2 = clat * clat
        clat4 = clat2 * clat2
        slat2 = slat * slat

        plg[0][0] = 1.0
        plg[1][0] = clat
        plg[2][0] = 0.5 * (3.0 * clat2 - 1.0)
        plg[3][0] = 0.5 * (5.0 * clat * clat2 - 3.0 * clat)
        plg[4][0] = (35.0 * clat4 - 30.0 * clat2 + 3.0) / 8.0
        plg[5][0] = (
            63.0 * clat2 * clat2 * clat
            - 70.0 * clat2 * clat
            + 15.0 * clat
        ) / 8.0
        plg[6][0] = (11.0 * clat * plg[5][0] - 5.0 * plg[4][0]) / 6.0

        plg[1][1] = slat
        plg[2][1] = 3.0 * clat * slat
        plg[3][1] = 1.5 * (5.0 * clat2 - 1.0) * slat
        plg[4][1] = 2.5 * (7.0 * clat2 * clat - 3.0 * clat) * slat
        plg[5][1] = 1.875 * (21.0 * clat4 - 14.0 * clat2 + 1.0) * slat
        plg[6][1] = (11.0 * clat * plg[5][1] - 6.0 * plg[4][1]) / 5.0

        plg[2][2] = 3.0 * slat2
        plg[3][2] = 15.0 * slat2 * clat
        plg[4][2] = 7.5 * (7.0 * clat2 - 1.0) * slat2
        plg[5][2] = 3.0 * clat * plg[4][2] - 2.0 * plg[3][2]
        plg[6][2] = (11.0 * clat * plg[5][2] - 7.0 * plg[4][2]) / 4.0

        plg[3][3] = 15.0 * slat2 * slat
        plg[4][3] = 105.0 * slat2 * slat * clat
        plg[5][3] = (9.0 * clat * plg[4][3] - 7.0 * plg[3][3]) / 2.0
        plg[6][3] = (11.0 * clat * plg[5][3] - 8.0 * plg[4][3]) / 3.0
        lastlat = lat

    if doy != lastdoy:
        cdoy[1] = math.cos(doy2rad * doy)
        sdoy[1] = math.sin(doy2rad * doy)
        cdoy[2] = math.cos(doy2rad * doy * 2.0)
        sdoy[2] = math.sin(doy2rad * doy * 2.0)
        lastdoy = doy

    lst = _fortran_mod(utsec / 3600.0 + lon / 15.0 + 24.0, 24.0)
    if lst != lastlst:
        clst[1] = math.cos(lst2rad * lst)
        slst[1] = math.sin(lst2rad * lst)
        clst[2] = math.cos(lst2rad * lst * 2.0)
        slst[2] = math.sin(lst2rad * lst * 2.0)
        clst[3] = math.cos(lst2rad * lst * 3.0)
        slst[3] = math.sin(lst2rad * lst * 3.0)
        lastlst = lst

    if lon != lastlon:
        clon[1] = math.cos(deg2rad * lon)
        slon[1] = math.sin(deg2rad * lon)
        clon[2] = math.cos(deg2rad * lon * 2.0)
        slon[2] = math.sin(deg2rad * lon * 2.0)
        lastlon = lon

    # Coupled linear terms.
    c = ctimeind
    for n in range(0, amaxn + 1):
        bf[c] = plg[n][0]
        c += 1

    if c != cintann:
        raise RuntimeError("problem with basis definitions")
    for s in range(1, amaxs + 1):
        cosdoy = cdoy[s]
        sindoy = sdoy[s]
        for n in range(0, amaxn + 1):
            pl = plg[n][0]
            bf[c] = pl * cosdoy
            bf[c + 1] = pl * sindoy
            c += 2

    if c != ctide:
        raise RuntimeError("problem with basis definitions")
    for l in range(1, tmaxl + 1):
        coslst = clst[l]
        sinlst = slst[l]
        for n in range(l, tmaxn + 1):
            pl = plg[n][l]
            bf[c] = pl * coslst
            bf[c + 1] = pl * sinlst
            c += 2
        for s in range(1, tmaxs + 1):
            cosdoy = cdoy[s]
            sindoy = sdoy[s]
            for n in range(l, tmaxn + 1):
                pl = plg[n][l]
                bf[c] = pl * coslst * cosdoy
                bf[c + 1] = pl * sinlst * cosdoy
                bf[c + 2] = pl * coslst * sindoy
                bf[c + 3] = pl * sinlst * sindoy
                c += 4

    if c != cspw:
        raise RuntimeError("problem with basis definitions")
    for m in range(1, pmaxm + 1):
        coslon = clon[m]
        sinlon = slon[m]
        for n in range(m, pmaxn + 1):
            pl = plg[n][m]
            bf[c] = pl * coslon
            bf[c + 1] = pl * sinlon
            c += 2
        for s in range(1, pmaxs + 1):
            cosdoy = cdoy[s]
            sindoy = sdoy[s]
            for n in range(m, pmaxn + 1):
                pl = plg[n][m]
                bf[c] = pl * coslon * cosdoy
                bf[c + 1] = pl * sinlon * cosdoy
                bf[c + 2] = pl * coslon * sindoy
                bf[c + 3] = pl * sinlon * sindoy
                c += 4

    if c != csfx:
        raise RuntimeError("problem with basis definitions")
    dfa = sfluxavg - sfluxavgref
    df = sflux - sfluxavg
    bf[c] = dfa
    bf[c + 1] = dfa * dfa
    bf[c + 2] = df
    bf[c + 3] = df * df
    bf[c + 4] = df * dfa
    c += nsfx

    if c != cextra:
        raise RuntimeError("problem with basis definitions")
    sza = solzen(doy, lst, lat, lon)
    bf[c] = -0.5 * math.tanh((sza - 98.0) / 6.0)
    bf[c + 1] = -0.5 * math.tanh((sza - 101.5) / 20.0)
    bf[c + 2] = dfa * bf[c]
    bf[c + 3] = dfa * bf[c + 1]
    bf[c + 4] = dfa * plg[2][0]
    bf[c + 5] = dfa * plg[4][0]
    bf[c + 6] = dfa * plg[0][0] * cdoy[1]
    bf[c + 7] = dfa * plg[0][0] * sdoy[1]
    bf[c + 8] = dfa * plg[0][0] * cdoy[2]
    bf[c + 9] = dfa * plg[0][0] * sdoy[2]
    if sfluxavg <= sfluxavg_quad_cutoff:
        bf[c + 10] = dfa * dfa
    else:
        cutoff_delta = sfluxavg_quad_cutoff - sfluxavgref
        bf[c + 10] = cutoff_delta * (2.0 * dfa - cutoff_delta)
    bf[c + 11] = bf[c + 10] * plg[2][0]
    bf[c + 12] = bf[c + 10] * plg[4][0]
    bf[c + 13] = df * plg[2][0]
    bf[c + 14] = df * plg[4][0]

    # Nonlinear terms.
    c = cnonlin
    if c != csfxmod:
        raise RuntimeError("problem with basis definitions")
    bf[c] = dfa
    bf[c + 1] = dfa * dfa
    bf[c + 2] = df
    bf[c + 3] = df * df
    bf[c + 4] = df * dfa
    c += nsfxmod

    if c != cmag:
        raise RuntimeError("problem with basis set")
    for j in range(7):
        bf[c + j] = float(ap[j]) - _f32(4.0)
    # c+7 intentionally remains zero.
    bf[c + 8] = doy2rad * doy
    bf[c + 9] = lst2rad * lst
    bf[c + 10] = deg2rad * lon
    bf[c + 11] = lst2rad * utsec / _f32(3600.0)
    bf[c + 12] = abs(lat)
    c += 13
    for m in range(0, 2):
        for n in range(0, amaxn + 1):
            bf[c] = plg[n][m]
            c += 1

    c = cut
    bf[c] = lst2rad * utsec / _f32(3600.0)
    bf[c + 1] = doy2rad * doy
    bf[c + 2] = dfa
    bf[c + 3] = deg2rad * lon
    bf[c + 4] = plg[1][0]
    bf[c + 5] = plg[3][0]
    bf[c + 6] = plg[5][0]
    bf[c + 7] = plg[3][2]
    bf[c + 8] = plg[5][2]

    swg = _parameters.swg
    for j in range(0, mbf + 1):
        if not swg[j]:
            bf[j] = 0.0

    return bf


def solzen(ddd: float, lst: float, lat: float, lon: float) -> float:
    """Calculate solar zenith angle using the authoritative source semantics."""
    # ``wlon`` and the first ``teqnx`` assignment are retained for source
    # traceability, even though the latter is immediately overwritten.
    wlon = _f32(360.0) - float(lon)
    teqnx = float(ddd) + (float(lst) + wlon / 15.0) / 24.0 + 0.9369
    del teqnx
    teqnx = float(ddd) + 0.9369

    p1, p2, p3, p4, p5 = _SOLZEN_P
    dec = (
        23.256 * math.sin(p1 * (teqnx - 82.242))
        + 0.381 * math.sin(p2 * (teqnx - 44.855))
        + 0.167 * math.sin(p3 * (teqnx - 23.355))
        - 0.013 * math.sin(p4 * (teqnx + 11.97))
        + 0.011 * math.sin(p5 * (teqnx - 10.410))
        + 0.339137
    )
    dec = dec * deg2rad

    tf = teqnx - 0.5
    teqt = (
        -7.38 * math.sin(p1 * (tf - 4.0))
        - 9.87 * math.sin(p2 * (tf + 9.0))
        + 0.27 * math.sin(p3 * (tf - 53.0))
        - 0.2 * math.cos(p4 * (tf - 17.0))
    )

    phi = (pi / 12.0) * (float(lst) - 12.0) + teqt * deg2rad / 4.0
    rlat = float(lat) * deg2rad
    cosx = (
        math.sin(rlat) * math.sin(dec)
        + math.cos(rlat) * math.cos(dec) * math.cos(phi)
    )
    if abs(cosx) > 1.0:
        cosx = math.copysign(1.0, cosx)
    return math.acos(cosx) / deg2rad


def sfluxmod(iz: int, gf: Sequence[float], parmset, dffact: float) -> float:
    """Legacy nonlinear solar-flux modulation of annual/tidal/SPW terms."""
    _require_len(gf, maxnbf, "gf")
    if parmset.beta is None:
        raise ValueError("parmset.beta is not allocated")

    swg = _parameters.swg
    if swg[csfxmod]:
        f1 = parmset.beta[csfxmod, iz] * gf[csfxmod] + (
            parmset.beta[csfx + 2, iz] * gf[csfxmod + 2]
            + parmset.beta[csfx + 3, iz] * gf[csfxmod + 3]
        ) * dffact
    else:
        f1 = 0.0

    if swg[csfxmod + 1]:
        f2 = parmset.beta[csfxmod + 1, iz] * gf[csfxmod] + (
            parmset.beta[csfx + 2, iz] * gf[csfxmod + 2]
            + parmset.beta[csfx + 3, iz] * gf[csfxmod + 3]
        ) * dffact
    else:
        f2 = 0.0

    if swg[csfxmod + 2]:
        f3 = parmset.beta[csfxmod + 2, iz] * gf[csfxmod]
    else:
        f3 = 0.0

    total = _f32(0.0)
    for j in range(0, mbf + 1):
        if _parameters.zsfx[j]:
            total = total + parmset.beta[j, iz] * gf[j] * f1
            continue
        if _parameters.tsfx[j]:
            total = total + parmset.beta[j, iz] * gf[j] * f2
            continue
        if _parameters.psfx[j]:
            total = total + parmset.beta[j, iz] * gf[j] * f3
            continue
    return total


def _g0fn(a: float, k00r: float, k00s: float) -> float:
    return a + (k00r - 1.0) * (a + (math.exp(-a * k00s) - 1.0) / k00s)


def geomag(p0: Sequence[float], bf: Sequence[float], plg_input: Sequence[Sequence[float]]) -> float:
    """Legacy nonlinear geomagnetic dependence."""
    _require_len(p0, nmag, "p0")
    _require_len(bf, 13, "bf")
    if len(plg_input) != 7 or any(len(row) != 2 for row in plg_input):
        raise ValueError("plg_input must have logical shape (0:6,0:1)")

    swg = _parameters.swg
    if not (swg[cmag] or swg[cmag + 1]):
        return 0.0

    p = [float(x) for x in p0]
    swg1 = list(swg[cmag : cmag + nmag])

    if swg1[0] == swg1[1]:
        if p[1] == 0.0:
            return 0.0
        for j in range(2, 26):
            if not swg1[j]:
                p[j] = 0.0
        p[8] = float(p0[8])
        del_a = _g0fn(float(bf[0]), p[0], p[1])
        value = (
            p[2] * plg_input[0][0]
            + p[3] * plg_input[2][0]
            + p[4] * plg_input[4][0]
            + (p[5] * plg_input[1][0] + p[6] * plg_input[3][0] + p[7] * plg_input[5][0])
            * math.cos(float(bf[8]) - p[8])
            + (p[9] * plg_input[1][1] + p[10] * plg_input[3][1] + p[11] * plg_input[5][1])
            * math.cos(float(bf[9]) - p[12])
            + (1.0 + p[13] * plg_input[1][0])
            * (p[14] * plg_input[2][1] + p[15] * plg_input[4][1] + p[16] * plg_input[6][1])
            * math.cos(float(bf[10]) - p[17])
            + (p[18] * plg_input[1][1] + p[19] * plg_input[3][1] + p[20] * plg_input[5][1])
            * math.cos(float(bf[10]) - p[21])
            * math.cos(float(bf[8]) - p[8])
            + (p[22] * plg_input[1][0] + p[23] * plg_input[3][0] + p[24] * plg_input[5][0])
            * math.cos(float(bf[11]) - p[25])
        ) * del_a
        return value

    if p[28] == 0.0:
        return _f32(0.0)
    for j in range(30, nmag):
        if not swg1[j]:
            p[j] = _f32(0.0)
    p[36] = float(p0[36])
    gbeta = p[28] / (1.0 + p[29] * (45.0 - float(bf[12])))
    ex = math.exp(-10800.0 * gbeta)
    sumex = 1.0 + (1.0 - ex**19.0) * ex**0.5 / (1.0 - ex)
    g = [0.0] * 7
    for i in range(1, 7):
        g[i] = _g0fn(float(bf[i]), p[26], p[27])
    del_a = (
        g[1]
        + (
            g[2] * ex
            + g[3] * ex * ex
            + g[4] * ex**3.0
            + (g[5] * ex**4.0 + g[6] * ex**12.0) * (1.0 - ex**8.0) / (1.0 - ex)
        )
    ) / sumex
    value = (
        p[30] * plg_input[0][0]
        + p[31] * plg_input[2][0]
        + p[32] * plg_input[4][0]
        + (p[33] * plg_input[1][0] + p[34] * plg_input[3][0] + p[35] * plg_input[5][0])
        * math.cos(float(bf[8]) - p[36])
        + (p[37] * plg_input[1][1] + p[38] * plg_input[3][1] + p[39] * plg_input[5][1])
        * math.cos(float(bf[9]) - p[40])
        + (1.0 + p[41] * plg_input[1][0])
        * (p[42] * plg_input[2][1] + p[43] * plg_input[4][1] + p[44] * plg_input[6][1])
        * math.cos(float(bf[10]) - p[45])
        + (p[46] * plg_input[1][1] + p[47] * plg_input[3][1] + p[48] * plg_input[5][1])
        * math.cos(float(bf[10]) - p[49])
        * math.cos(float(bf[8]) - p[36])
        + (p[50] * plg_input[1][0] + p[51] * plg_input[3][0] + p[52] * plg_input[5][0])
        * math.cos(float(bf[11]) - p[53])
    ) * del_a
    return value


def utdep(p0: Sequence[float], bf: Sequence[float]) -> float:
    """Legacy nonlinear universal-time dependence."""
    _require_len(p0, nut, "p0")
    _require_len(bf, 9, "bf")
    p = [float(x) for x in p0]
    swg1 = list(_parameters.swg[cut : cut + nut])
    for j in range(3, nut):
        if not swg1[j]:
            p[j] = _f32(0.0)
    return (
        math.cos(float(bf[0]) - p[0])
        * (1.0 + p[3] * float(bf[4]) * math.cos(float(bf[1]) - p[1]))
        * (1.0 + p[4] * float(bf[2]))
        * (1.0 + p[5] * float(bf[4]))
        * (p[6] * float(bf[4]) + p[7] * float(bf[5]) + p[8] * float(bf[6]))
        + math.cos(float(bf[0]) - p[2] + 2 * float(bf[3]))
        * (p[9] * float(bf[7]) + p[10] * float(bf[8]))
        * (1.0 + p[11] * float(bf[2]))
    )


__all__ = [
    "plg",
    "cdoy",
    "sdoy",
    "clst",
    "slst",
    "clon",
    "slon",
    "sfluxavgref",
    "sfluxavg_quad_cutoff",
    "lastlat",
    "lastdoy",
    "lastlst",
    "lastlon",
    "globe",
    "solzen",
    "sfluxmod",
    "geomag",
    "utdep",
]
