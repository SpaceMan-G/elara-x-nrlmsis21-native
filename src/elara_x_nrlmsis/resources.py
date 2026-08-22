"""Verified external parameter-resource handling for native NRLMSIS 2.1.

This module does not redistribute the official ``msis21.parm`` payload.
Resolution precedence:
1. explicit resource file supplied by the caller;
2. ``ELARA_X_NRLMSIS21_PARM`` environment variable.

There is no implicit current-working-directory lookup and no network
acquisition in this resolver.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
from typing import Mapping, Optional, Union

RESOURCE_BASENAME = "msis21.parm"
RESOURCE_SHA256 = "a322a749f368e73117dd20f3fdcf7389dabc5509f4c27073cc5580999381b508"
RESOURCE_BYTES = 536576
RESOURCE_SHAPE = (512, 131)
RESOURCE_SCALAR_COUNT = 67072
RESOURCE_ENDIANNESS = "little"
RESOURCE_ENVVAR = "ELARA_X_NRLMSIS21_PARM"

PathLike = Union[str, os.PathLike[str]]


class ResourceError(RuntimeError):
    """Base class for controlled NRLMSIS 2.1 resource failures."""


class ResourceNotConfiguredError(ResourceError):
    """No explicit resource and no configured environment resource."""


class ResourceNotFoundError(ResourceError):
    """The configured resource path does not identify an existing regular file."""


class ResourceIdentityError(ResourceError):
    """The configured file does not match the frozen official resource identity."""


class ResourceInitializationError(ResourceError):
    """The verified resource could not be loaded by the frozen scientific loader."""


@dataclass(frozen=True)
class VerifiedParameterResource:
    path: Path
    sha256: str
    bytes: int
    basename: str = RESOURCE_BASENAME
    shape: tuple[int, int] = RESOURCE_SHAPE
    scalar_count: int = RESOURCE_SCALAR_COUNT
    endianness: str = RESOURCE_ENDIANNESS


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_parameter_resource(
    resource_file: Optional[PathLike] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> Path:
    """Resolve only the two routes frozen by the Stage-9 contract."""
    if resource_file is not None:
        text = os.fspath(resource_file)
        if not str(text).strip():
            raise ResourceNotConfiguredError("Explicit NRLMSIS 2.1 resource path is empty.")
        return Path(text).expanduser().resolve(strict=False)

    env = os.environ if environ is None else environ
    configured = env.get(RESOURCE_ENVVAR)
    if configured is not None and str(configured).strip():
        return Path(configured).expanduser().resolve(strict=False)

    raise ResourceNotConfiguredError(
        "NRLMSIS 2.1 parameter resource is not configured. Supply resource_file="
        " explicitly or set ELARA_X_NRLMSIS21_PARM."
    )


def verify_parameter_resource(resource_file: PathLike) -> VerifiedParameterResource:
    """Verify exact official resource identity before scientific loading."""
    path = Path(resource_file).expanduser().resolve(strict=False)

    if path.name != RESOURCE_BASENAME:
        raise ResourceIdentityError(
            f"Expected resource basename {RESOURCE_BASENAME!r}; got {path.name!r}."
        )
    if not path.is_file():
        raise ResourceNotFoundError(f"NRLMSIS 2.1 parameter resource not found: {path}")

    size = path.stat().st_size
    if size != RESOURCE_BYTES:
        raise ResourceIdentityError(
            f"NRLMSIS 2.1 resource byte count mismatch: expected {RESOURCE_BYTES}, got {size}."
        )

    digest = _sha256_file(path)
    if digest != RESOURCE_SHA256:
        raise ResourceIdentityError(
            "NRLMSIS 2.1 resource SHA-256 mismatch: "
            f"expected {RESOURCE_SHA256}, got {digest}."
        )

    return VerifiedParameterResource(path=path, sha256=digest, bytes=size)


def resolve_and_verify_parameter_resource(
    resource_file: Optional[PathLike] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> VerifiedParameterResource:
    return verify_parameter_resource(
        resolve_parameter_resource(resource_file, environ=environ)
    )


def initialize_nrlmsis21(
    resource_file: Optional[PathLike] = None,
    *,
    environ: Optional[Mapping[str, str]] = None,
) -> VerifiedParameterResource:
    """Verify, load through frozen ``parameters.msisinit``, then reverify."""
    verified_before = resolve_and_verify_parameter_resource(
        resource_file, environ=environ
    )

    from . import parameters

    try:
        parameters.msisinit(
            parmpath=str(verified_before.path.parent) + os.sep,
            parmfile=verified_before.path.name,
        )
    except Exception as exc:
        raise ResourceInitializationError(
            f"Frozen NRLMSIS 2.1 loader failed for verified resource: {verified_before.path}"
        ) from exc

    verified_after = verify_parameter_resource(verified_before.path)
    if (
        verified_after.sha256 != verified_before.sha256
        or verified_after.bytes != verified_before.bytes
    ):
        raise ResourceIdentityError(
            "NRLMSIS 2.1 parameter resource changed during initialization."
        )
    return verified_after
