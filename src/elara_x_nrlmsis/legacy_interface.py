"""
Native NRLMSIS 2.1 legacy GTD8D compatibility interface.

Authoritative counterpart
-------------------------
NRL NRLMSIS 2.1 ``msis_gtd8d.F90``.

Derivative translation notice
-----------------------------
This file translates the historical GTD8D wrapper semantics around the native
NRLMSIS 2.1 ``msiscalc`` implementation.  It intentionally preserves the
legacy REAL(4) input/output boundary, Fortran MOD semantics for YYDDD input,
ignored STL/MASS arguments, SI-to-CGS density conversion, missing-value guard,
and the historical species ordering.

Parameter-resource acquisition/distribution remains a separate Stage-9 concern;
this wrapper delegates initialization behavior to ``model.msiscalc`` exactly as
the authoritative GTD8D wrapper delegates to MSISCALC.
"""

from __future__ import annotations

import math
import struct
from typing import Sequence

from .constants import dmissing
from .model import msiscalc


def _f32(value: float) -> float:
    """Round a value to IEEE binary32 and promote the result to Python float."""
    return struct.unpack("<f", struct.pack("<f", float(value)))[0]


def _fortran_mod_int(a: int, p: int) -> int:
    """Fortran MOD for integer operands (quotient truncated toward zero)."""
    a = int(a)
    p = int(p)
    if p == 0:
        raise ZeroDivisionError("Fortran MOD divisor must be non-zero")
    qmag = abs(a) // abs(p)
    q = -qmag if (a < 0) != (p < 0) else qmag
    return a - q * p


def _legacy_ap(ap: Sequence[float]) -> tuple[float, ...]:
    if len(ap) != 7:
        raise ValueError("ap must contain exactly 7 values")
    return tuple(_f32(x) for x in ap)


def gtd8d(
    iyd: int,
    sec: float,
    alt: float,
    glat: float,
    glong: float,
    stl: float,
    f107a: float,
    f107: float,
    ap: Sequence[float],
    mass: int,
):
    """Evaluate the NRLMSIS 2.1 GTD8D legacy compatibility interface.

    Inputs corresponding to Fortran ``REAL(4)`` are first rounded to binary32
    before promotion to the native model's binary64 arithmetic.  ``stl`` and
    ``mass`` are intentionally accepted but ignored, matching the authoritative
    wrapper.  The return value is ``(d, t)`` where ``d`` is the ten-element
    legacy density tuple and ``t=(T_exo, T_alt)``; every returned value has been
    rounded to the authoritative REAL(4) output boundary.
    """

    # Authoritative REAL(4) argument boundary.  STL is deliberately ignored
    # after arrival at that boundary, and MASS is not forwarded by GTD8D.
    xday = float(_fortran_mod_int(int(iyd), 1000))
    xutsec = _f32(sec)
    xalt = _f32(alt)
    xlat = _f32(glat)
    xlon = _f32(glong)
    _ = _f32(stl)
    xsfluxavg = _f32(f107a)
    xsflux = _f32(f107)
    xap = _legacy_ap(ap)
    _ = int(mass)

    xtn, xdn_tuple, xtex = msiscalc(
        xday,
        xutsec,
        xalt,
        xlat,
        xlon,
        xsfluxavg,
        xsflux,
        xap,
        return_tex=True,
    )

    # T(1)=TEX and T(2)=local temperature, both rounded through SNGL.
    t = (_f32(xtex), _f32(xtn))

    # Authoritative WHERE guard: the missing sentinel is not unit-scaled.
    xdn = [float(v) for v in xdn_tuple]
    for i, value in enumerate(xdn):
        if value != dmissing:
            xdn[i] = value * 1.0e-6
    # Total mass density has a further x1000 conversion, yielding net 1e-3
    # from kg/m^3 to g/cm^3 after the common 1e-6 operation.
    if xdn[0] != dmissing:
        xdn[0] = xdn[0] * 1.0e3

    # Exact legacy species ordering followed by the REAL(4)/SNGL boundary.
    d = (
        _f32(xdn[4]),  # He
        _f32(xdn[3]),  # O
        _f32(xdn[1]),  # N2
        _f32(xdn[2]),  # O2
        _f32(xdn[6]),  # Ar
        _f32(xdn[0]),  # total mass density
        _f32(xdn[5]),  # H
        _f32(xdn[7]),  # N
        _f32(xdn[8]),  # anomalous O
        _f32(xdn[9]),  # NO
    )
    return d, t


__all__ = ["gtd8d"]
