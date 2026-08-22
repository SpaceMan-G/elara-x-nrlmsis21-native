from __future__ import annotations

import hashlib
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import elara_x_nrlmsis.resources as r


class TestVerifiedExternalResourceResolver(unittest.TestCase):
    def _fake_resource(self, directory: Path, payload: bytes = b"abc") -> Path:
        path = directory / r.RESOURCE_BASENAME
        path.write_bytes(payload)
        return path

    def test_01_resource_identity_constants_are_frozen(self):
        self.assertEqual(r.RESOURCE_BASENAME, "msis21.parm")
        self.assertEqual(r.RESOURCE_SHA256, "a322a749f368e73117dd20f3fdcf7389dabc5509f4c27073cc5580999381b508")
        self.assertEqual(r.RESOURCE_BYTES, 536576)
        self.assertEqual(r.RESOURCE_SHAPE, (512, 131))
        self.assertEqual(r.RESOURCE_SCALAR_COUNT, 67072)
        self.assertEqual(r.RESOURCE_ENDIANNESS, "little")

    def test_02_explicit_path_has_precedence_over_environment(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            explicit = root / "explicit" / "msis21.parm"
            envpath = root / "env" / "msis21.parm"
            got = r.resolve_parameter_resource(
                explicit,
                environ={r.RESOURCE_ENVVAR: str(envpath)},
            )
            self.assertEqual(got, explicit.resolve())

    def test_03_environment_path_is_second_resolution_route(self):
        with tempfile.TemporaryDirectory() as td:
            envpath = Path(td) / "env" / "msis21.parm"
            got = r.resolve_parameter_resource(
                environ={r.RESOURCE_ENVVAR: str(envpath)}
            )
            self.assertEqual(got, envpath.resolve())

    def test_04_no_configuration_fails_closed(self):
        with self.assertRaises(r.ResourceNotConfiguredError):
            r.resolve_parameter_resource(environ={})

    def test_05_empty_explicit_path_does_not_become_cwd(self):
        with self.assertRaises(r.ResourceNotConfiguredError):
            r.resolve_parameter_resource("", environ={})

    def test_06_no_implicit_cwd_fallback(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / r.RESOURCE_BASENAME).write_bytes(b"x")
            old = Path.cwd()
            try:
                os.chdir(root)
                with self.assertRaises(r.ResourceNotConfiguredError):
                    r.resolve_parameter_resource(environ={})
            finally:
                os.chdir(old)

    def test_07_wrong_basename_fails_before_loading(self):
        with tempfile.TemporaryDirectory() as td:
            bad = Path(td) / "renamed.parm"
            bad.write_bytes(b"x")
            with self.assertRaises(r.ResourceIdentityError):
                r.verify_parameter_resource(bad)

    def test_08_missing_resource_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(r.ResourceNotFoundError):
                r.verify_parameter_resource(Path(td) / r.RESOURCE_BASENAME)

    def test_09_wrong_byte_count_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / r.RESOURCE_BASENAME
            path.write_bytes(b"x")
            with self.assertRaisesRegex(r.ResourceIdentityError, "byte count"):
                r.verify_parameter_resource(path)

    def test_10_wrong_sha256_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / r.RESOURCE_BASENAME
            path.write_bytes(b"\x00" * r.RESOURCE_BYTES)
            with self.assertRaisesRegex(r.ResourceIdentityError, "SHA-256"):
                r.verify_parameter_resource(path)

    def test_11_verified_resource_metadata(self):
        with tempfile.TemporaryDirectory() as td:
            payload = b"verified"
            path = self._fake_resource(Path(td), payload)
            with mock.patch.object(r, "RESOURCE_BYTES", len(payload)),                  mock.patch.object(r, "RESOURCE_SHA256", hashlib.sha256(payload).hexdigest()):
                verified = r.verify_parameter_resource(path)
            self.assertEqual(verified.path, path.resolve())
            self.assertEqual(verified.bytes, len(payload))
            self.assertEqual(verified.sha256, hashlib.sha256(payload).hexdigest())
            self.assertEqual(verified.basename, "msis21.parm")

    def test_12_resolve_and_verify_uses_environment(self):
        with tempfile.TemporaryDirectory() as td:
            payload = b"environment"
            path = self._fake_resource(Path(td), payload)
            with mock.patch.object(r, "RESOURCE_BYTES", len(payload)),                  mock.patch.object(r, "RESOURCE_SHA256", hashlib.sha256(payload).hexdigest()):
                verified = r.resolve_and_verify_parameter_resource(
                    environ={r.RESOURCE_ENVVAR: str(path)}
                )
            self.assertEqual(verified.path, path.resolve())

    def test_13_initialize_calls_frozen_loader_with_parent_separator(self):
        with tempfile.TemporaryDirectory() as td:
            payload = b"loader"
            path = self._fake_resource(Path(td), payload)
            from elara_x_nrlmsis import parameters
            called = {}

            def fake_msisinit(*, parmpath="", parmfile="msis21.parm", **kwargs):
                called["parmpath"] = parmpath
                called["parmfile"] = parmfile

            with mock.patch.object(r, "RESOURCE_BYTES", len(payload)),                  mock.patch.object(r, "RESOURCE_SHA256", hashlib.sha256(payload).hexdigest()),                  mock.patch.object(parameters, "msisinit", fake_msisinit):
                verified = r.initialize_nrlmsis21(path)

            self.assertEqual(called["parmpath"], str(path.resolve().parent) + os.sep)
            self.assertEqual(called["parmfile"], "msis21.parm")
            self.assertEqual(verified.path, path.resolve())

    def test_14_initialize_wraps_loader_failure(self):
        with tempfile.TemporaryDirectory() as td:
            payload = b"failure"
            path = self._fake_resource(Path(td), payload)
            from elara_x_nrlmsis import parameters

            def fail(**kwargs):
                raise RuntimeError("loader failure")

            with mock.patch.object(r, "RESOURCE_BYTES", len(payload)),                  mock.patch.object(r, "RESOURCE_SHA256", hashlib.sha256(payload).hexdigest()),                  mock.patch.object(parameters, "msisinit", fail):
                with self.assertRaises(r.ResourceInitializationError):
                    r.initialize_nrlmsis21(path)

    def test_15_initialize_reverifies_resource_after_load(self):
        with tempfile.TemporaryDirectory() as td:
            payload = b"stable"
            path = self._fake_resource(Path(td), payload)
            from elara_x_nrlmsis import parameters

            def mutate(**kwargs):
                path.write_bytes(b"changed")

            with mock.patch.object(r, "RESOURCE_BYTES", len(payload)),                  mock.patch.object(r, "RESOURCE_SHA256", hashlib.sha256(payload).hexdigest()),                  mock.patch.object(parameters, "msisinit", mutate):
                with self.assertRaises(r.ResourceIdentityError):
                    r.initialize_nrlmsis21(path)


if __name__ == "__main__":
    unittest.main()
