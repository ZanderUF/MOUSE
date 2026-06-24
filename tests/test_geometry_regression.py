"""Characterization tests pinning drum/vessel/reflector geometry outputs.

These lock the numeric results of the reactivity-control and vessel-stack
calculations for the three production reactors so the capability-driven refactor
(control-element strategy + vessel-stack model) stays byte-for-byte
behaviour-preserving. The golden values were captured from the pre-refactor
implementation.
"""

from __future__ import annotations

import unittest

import reactor_test_support as support


# Golden values captured from the pre-refactor code for the canonical design
# points in reactor_test_support.REACTOR_BUILD_CASES.
GOLDEN = {
    "LTMR": {
        "Drum Count": 12,
        "Drum Radius": 9.39313130825758,
        "All Drums Area": 3326.22716178754,
        "Control Drum Absorber Mass": 76.08921100882458,
        "Control Drum Reflector Mass": 670.6413363424225,
        "Control Drums Mass": 746.7305473512471,
        "Core Radius": 52.94604756474369,
        "Radial Reflector Mass": 374.12805148854414,
        "Axial Reflector Mass": 282.46139068726467,
        "Moderator Mass": 382.83481161247545,
        "Vessel Mass": 1775.2861362983704,
        "Guard Vessel Mass": 1002.9391573624496,
        "Cooling Vessel Mass": 545.8938785995205,
        "Intake Vessel Mass": 592.6288397323277,
        "Vessels Total Radius": 80.60604756474369,
        "Vessels Total Height": 282.8724849499458,
        "Total Vessels Mass": 3916.7480119926686,
    },
    "GCMR": {
        "Drum Count": 24,
        "Drum Radius": 9.530986101432001,
        "All Drums Area": 6849.15172354323,
        "Control Drum Absorber Mass": 241.92809928478314,
        "Control Drum Reflector Mass": 2163.98968500009,
        "Control Drums Mass": 2405.917784284873,
        "Core Radius": 107.17064371832429,
        "Radial Reflector Mass": 2818.80505680517,
        "Axial Reflector Mass": 1124.9537914668945,
        "Moderator Mass": 3983.3799067517043,
        "Vessel Mass": 6775.0160187052725,
        "Guard Vessel Mass": 0.0,
        "Cooling Vessel Mass": 1041.5627994272832,
        "Intake Vessel Mass": 1101.4508801031109,
        "Vessels Total Radius": 120.67064371832429,
        "Vessels Total Height": 334.2665715851499,
        "Total Vessels Mass": 8918.029698235667,
        "RPV Outer Radius": 110.17064371832429,
        "RPV Outer Height": 358.89557158514987,
    },
    "HPMR": {
        "Drum Count": 12,
        "Drum Radius": 49.71794080574132,
        "All Drums Area": 93187.4407395295,
        "Control Drum Absorber Mass": 1047.7477663063564,
        "Control Drum Reflector Mass": 49450.055978822114,
        "Control Drums Mass": 50497.80374512847,
        "Core Radius": 243.9477060470781,
        "Radial Reflector Mass": 22506.69413677664,
        "Axial Reflector Mass": 67938.59546906262,
        "Moderator Mass": 6788.122589570194,
        "Vessel Mass": 13086.154363276377,
        "Guard Vessel Mass": 0.0,
        "Cooling Vessel Mass": 3392.4754174790055,
        "Intake Vessel Mass": 3487.1869248319763,
        "Vessels Total Radius": 256.4477060470781,
        "Vessels Total Height": 447.8998767642017,
        "Total Vessels Mass": 19965.816705587356,
    },
}


class GeometryRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        support.install()
        cls.params = {rt: support.build_case(rt) for rt in GOLDEN}

    def test_golden_geometry_outputs_unchanged(self):
        for reactor_type, expected in GOLDEN.items():
            params = self.params[reactor_type]
            for key, value in expected.items():
                with self.subTest(reactor=reactor_type, key=key):
                    self.assertAlmostEqual(
                        float(params[key]), float(value), places=6
                    )

    def test_drum_reactors_publish_generic_control_element_keys(self):
        for reactor_type in GOLDEN:
            params = self.params[reactor_type]
            with self.subTest(reactor=reactor_type):
                self.assertEqual(params["Control Element Type"], "drum")
                self.assertEqual(params["Control Element Count"],
                                 params["Drum Count"])
                self.assertAlmostEqual(params["Control Element Mass"],
                                       params["Control Drums Mass"])

    def test_vessel_stack_breakdown_present(self):
        for reactor_type in GOLDEN:
            params = self.params[reactor_type]
            with self.subTest(reactor=reactor_type):
                stack = params["Vessel Stack"]
                self.assertEqual(len(stack), 4)  # concentric RVACS
                self.assertEqual(
                    sum(layer["mass"] for layer in stack),
                    params["Total Vessels Mass"],
                )


class BuildParamsDispatchTest(unittest.TestCase):
    """build_params distinguishes catalogued-but-unbuilt from truly unknown."""

    def setUp(self):
        support.install()
        from webapp.reactor_config import build_params
        self.build_params = build_params

    def test_registered_but_unimplemented_raises_not_implemented(self):
        # SFR is catalogued in the registry but has no params builder yet.
        with self.assertRaises(NotImplementedError):
            self.build_params("SFR", 10, 0.15, {}, active_height=100)

    def test_truly_unknown_type_raises_value_error(self):
        with self.assertRaises(ValueError):
            self.build_params("PWR", 10, 0.15, {})


if __name__ == "__main__":
    unittest.main()
