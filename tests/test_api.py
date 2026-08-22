from __future__ import annotations

import unittest
from unittest import mock

import elara_x_nrlmsis as ex
import elara_x_nrlmsis.api as api


EXPECTED_ALL = [
    "ModelNotInitializedError",
    "ResourceError",
    "ResourceNotConfiguredError",
    "ResourceNotFoundError",
    "ResourceIdentityError",
    "ResourceInitializationError",
    "VerifiedParameterResource",
    "initialize",
    "is_initialized",
    "calculate",
    "gtd8d",
]


class TestAssembledPublicAPI(unittest.TestCase):
    def test_01_package_all_exact(self):
        self.assertEqual(ex.__all__, EXPECTED_ALL)

    def test_02_package_root_exports_present(self):
        for name in EXPECTED_ALL:
            self.assertTrue(hasattr(ex, name), name)

    def test_03_raw_msiscalc_not_promoted_to_package_root(self):
        self.assertFalse(hasattr(ex, "msiscalc"))

    def test_04_is_initialized_false(self):
        with mock.patch.object(api._parameters, "initflag", False):
            self.assertFalse(api.is_initialized())

    def test_05_is_initialized_true(self):
        with mock.patch.object(api._parameters, "initflag", True):
            self.assertTrue(api.is_initialized())

    def test_06_calculate_fails_closed_before_initialization(self):
        with mock.patch.object(api._parameters, "initflag", False):
            with self.assertRaises(api.ModelNotInitializedError):
                api.calculate(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, [1]*7)

    def test_07_gtd8d_fails_closed_before_initialization(self):
        with mock.patch.object(api._parameters, "initflag", False):
            with self.assertRaises(api.ModelNotInitializedError):
                api.gtd8d(1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, [1]*7, 48)

    def test_08_initialize_delegates_explicit_resource(self):
        sentinel = object()
        with mock.patch.object(api._resources, "initialize_nrlmsis21", return_value=sentinel) as call:
            got = api.initialize("/tmp/msis21.parm")
        self.assertIs(got, sentinel)
        call.assert_called_once_with("/tmp/msis21.parm", environ=None)

    def test_09_initialize_delegates_environment_mapping(self):
        sentinel = object()
        env = {"ELARA_X_NRLMSIS21_PARM": "/tmp/msis21.parm"}
        with mock.patch.object(api._resources, "initialize_nrlmsis21", return_value=sentinel) as call:
            got = api.initialize(environ=env)
        self.assertIs(got, sentinel)
        call.assert_called_once_with(None, environ=env)

    def test_10_calculate_forwards_default_return_tex(self):
        sentinel = object()
        args = (172.0, 43200.0, 400.0, 45.0, 10.0, 150.0, 155.0, [4]*7)
        with mock.patch.object(api._parameters, "initflag", True),              mock.patch.object(api._model, "msiscalc", return_value=sentinel) as call:
            got = api.calculate(*args)
        self.assertIs(got, sentinel)
        call.assert_called_once_with(*args, return_tex=False)

    def test_11_calculate_forwards_return_tex_true(self):
        sentinel = object()
        args = (172.0, 43200.0, 400.0, 45.0, 10.0, 150.0, 155.0, [4]*7)
        with mock.patch.object(api._parameters, "initflag", True),              mock.patch.object(api._model, "msiscalc", return_value=sentinel) as call:
            got = api.calculate(*args, return_tex=True)
        self.assertIs(got, sentinel)
        call.assert_called_once_with(*args, return_tex=True)

    def test_12_gtd8d_forwards_all_arguments(self):
        sentinel = object()
        args = (24172, 43200.0, 400.0, 45.0, 10.0, 12.0, 150.0, 155.0, [4]*7, 48)
        with mock.patch.object(api._parameters, "initflag", True),              mock.patch.object(api._legacy_interface, "gtd8d", return_value=sentinel) as call:
            got = api.gtd8d(*args)
        self.assertIs(got, sentinel)
        call.assert_called_once_with(*args)

    def test_13_resource_error_propagates_from_initialize(self):
        err = api._resources.ResourceIdentityError("bad")
        with mock.patch.object(api._resources, "initialize_nrlmsis21", side_effect=err):
            with self.assertRaises(api._resources.ResourceIdentityError):
                api.initialize("/tmp/msis21.parm")

    def test_14_model_error_propagates_from_calculate(self):
        with mock.patch.object(api._parameters, "initflag", True),              mock.patch.object(api._model, "msiscalc", side_effect=ValueError("model")):
            with self.assertRaises(ValueError):
                api.calculate(1.0,2.0,3.0,4.0,5.0,6.0,7.0,[1]*7)

    def test_15_legacy_error_propagates_from_gtd8d(self):
        with mock.patch.object(api._parameters, "initflag", True),              mock.patch.object(api._legacy_interface, "gtd8d", side_effect=ValueError("legacy")):
            with self.assertRaises(ValueError):
                api.gtd8d(1,2.0,3.0,4.0,5.0,6.0,7.0,8.0,[1]*7,48)

    def test_16_repeated_initialize_has_no_api_local_cache(self):
        with mock.patch.object(api._resources, "initialize_nrlmsis21", return_value=object()) as call:
            api.initialize("/tmp/msis21.parm")
            api.initialize("/tmp/msis21.parm")
        self.assertEqual(call.call_count, 2)

    def test_17_model_not_initialized_error_is_runtime_error(self):
        self.assertTrue(issubclass(api.ModelNotInitializedError, RuntimeError))

    def test_18_model_not_initialized_message_points_to_initialize(self):
        with mock.patch.object(api._parameters, "initflag", False):
            try:
                api.calculate(1.0,2.0,3.0,4.0,5.0,6.0,7.0,[1]*7)
            except api.ModelNotInitializedError as exc:
                self.assertIn("initialize", str(exc))
            else:
                self.fail("ModelNotInitializedError was not raised")


if __name__ == "__main__":
    unittest.main()
