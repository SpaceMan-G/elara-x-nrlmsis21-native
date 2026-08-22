"""High-level, resource-safe public API for native NRLMSIS 2.1."""

from __future__ import annotations

from typing import Mapping, Optional, Sequence

from . import legacy_interface as _legacy_interface
from . import model as _model
from . import parameters as _parameters
from . import resources as _resources

PathLike = _resources.PathLike
VerifiedParameterResource = _resources.VerifiedParameterResource


class ModelNotInitializedError(RuntimeError):
    """Raised when assembled calculation is attempted before initialization."""


def is_initialized() -> bool:
    """Return whether the frozen native parameter state is initialized."""
    return bool(_parameters.initflag)


def initialize(
    resource_file: Optional[PathLike] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> VerifiedParameterResource:
    """Initialize through the Stage-9 verified external resource boundary."""
    return _resources.initialize_nrlmsis21(resource_file, environ=environ)


def _require_initialized() -> None:
    if not is_initialized():
        raise ModelNotInitializedError(
            "NRLMSIS 2.1 is not initialized. Call elara_x_nrlmsis.initialize(...) "
            "with the verified external msis21.parm resource before calculation."
        )


def calculate(
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
    """Run the frozen native model after verified initialization."""
    _require_initialized()
    return _model.msiscalc(
        day, utsec, z, lat, lon, sfluxavg, sflux, ap, return_tex=return_tex
    )


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
    """Run the frozen legacy GTD8D interface after verified initialization."""
    _require_initialized()
    return _legacy_interface.gtd8d(
        iyd, sec, alt, glat, glong, stl, f107a, f107, ap, mass
    )
