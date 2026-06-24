"""Tests for the control-element abstraction (drums vs. rods)."""

from __future__ import annotations

import math
import unittest

import reactor_test_support as support

support.install()  # stub openmc/watts before importing core_design modules

from core_design import control_elements as ce
from core_design.control_elements import (
    RodElement,
    DrumElement,
    resolve_control_element,
    _resolve_kind,
)


class DispatchTest(unittest.TestCase):
    def test_registry_drives_dispatch(self):
        self.assertEqual(_resolve_kind({"reactor type": "LTMR"}), "drum")
        self.assertEqual(_resolve_kind({"reactor type": "GCMR"}), "drum")
        self.assertEqual(_resolve_kind({"reactor type": "HPMR"}), "drum")
        self.assertEqual(_resolve_kind({"reactor type": "SFR"}), "rod")
        self.assertEqual(_resolve_kind({"reactor type": "SRE"}), "rod")

    def test_explicit_type_overrides_registry(self):
        self.assertEqual(
            _resolve_kind({"reactor type": "SFR", "Control Element Type": "drum"}),
            "drum",
        )
        self.assertEqual(_resolve_kind({"Control Element Type": "rod"}), "rod")

    def test_unknown_or_missing_reactor_falls_back_to_drum(self):
        self.assertEqual(_resolve_kind({}), "drum")
        self.assertEqual(_resolve_kind({"reactor type": "PWR"}), "drum")


class _DensityPatch:
    """Context manager that swaps collect_materials_data for a fixed table."""

    def __init__(self, densities):
        self._densities = {n: support.ThinMaterial(n, d)
                           for n, d in densities.items()}
        self._orig = None

    def __enter__(self):
        self._orig = ce.collect_materials_data
        ce.collect_materials_data = lambda params: self._densities
        return self

    def __exit__(self, *exc):
        ce.collect_materials_data = self._orig


class RodElementMathTest(unittest.TestCase):
    DENSITIES = {"B4C_enriched": 2.52, "SS304": 7.93}

    def _base_params(self, **overrides):
        p = {
            "Number of Control Rods": 7,
            "Control Rod Radius": 2.0,
            "Control Rod Clad Thickness": 0.2,
            "Control Rod Absorber": "B4C_enriched",
            "Control Rod Clad": "SS304",
            "Active Height": 100.0,
            "Axial Reflector Thickness": 20.0,
            "Control Rod Drive Mass": 5.0,
        }
        p.update(overrides)
        return p

    def test_masses_use_correct_radii_and_lengths(self):
        params = self._base_params()
        with _DensityPatch(self.DENSITIES):
            RodElement().compute(params)

        n, r_out, clad_t = 7, 2.0, 0.2
        r_abs = r_out - clad_t
        absorber_len = 100.0
        rod_len = 100.0 + 2 * 20.0  # default = absorber + 2*axial reflector

        exp_absorber = math.pi * r_abs**2 * absorber_len * n * 2.52 / 1000.0
        exp_clad = math.pi * (r_out**2 - r_abs**2) * rod_len * n * 7.93 / 1000.0
        exp_drive = n * 5.0

        self.assertAlmostEqual(params["Control Rod Absorber Mass"], exp_absorber)
        self.assertAlmostEqual(params["Control Rod Clad Mass"], exp_clad)
        self.assertAlmostEqual(params["Control Rod Drive Mass Total"], exp_drive)
        self.assertAlmostEqual(params["Control Rods Mass"], exp_absorber + exp_clad)
        self.assertAlmostEqual(
            params["Control Element Mass"], exp_absorber + exp_clad + exp_drive
        )

    def test_generic_keys_published(self):
        params = self._base_params()
        with _DensityPatch(self.DENSITIES):
            RodElement().compute(params)
        self.assertEqual(params["Control Element Type"], "rod")
        self.assertEqual(params["Control Element Count"], 7)
        self.assertAlmostEqual(params["All Rods Area"], math.pi * 2.0**2 * 7)
        self.assertAlmostEqual(
            params["All Rods Volume"], math.pi * 2.0**2 * (100 + 40) * 7
        )

    def test_explicit_absorber_and_rod_length_override_defaults(self):
        params = self._base_params(
            **{"Control Rod Absorber Length": 80.0, "Control Rod Length": 90.0}
        )
        with _DensityPatch(self.DENSITIES):
            RodElement().compute(params)
        r_abs = 1.8
        exp_absorber = math.pi * r_abs**2 * 80.0 * 7 * 2.52 / 1000.0
        self.assertAlmostEqual(params["Control Rod Absorber Mass"], exp_absorber)

    def test_zero_drive_mass_default(self):
        params = self._base_params()
        params.pop("Control Rod Drive Mass")
        with _DensityPatch(self.DENSITIES):
            RodElement().compute(params)
        self.assertEqual(params["Control Rod Drive Mass Total"], 0.0)
        self.assertAlmostEqual(
            params["Control Element Mass"], params["Control Rods Mass"]
        )


class RodElementValidationTest(unittest.TestCase):
    DENSITIES = {"B4C_enriched": 2.52, "SS304": 7.93}

    def _params(self, **overrides):
        p = {
            "Number of Control Rods": 7,
            "Control Rod Radius": 2.0,
            "Control Rod Clad Thickness": 0.2,
            "Control Rod Absorber": "B4C_enriched",
            "Control Rod Clad": "SS304",
            "Active Height": 100.0,
        }
        p.update(overrides)
        return p

    def test_clad_thicker_than_radius_rejected(self):
        with _DensityPatch(self.DENSITIES), self.assertRaises(ValueError):
            RodElement().compute(self._params(**{"Control Rod Clad Thickness": 2.5}))

    def test_nonpositive_rod_count_rejected(self):
        with _DensityPatch(self.DENSITIES), self.assertRaises(ValueError):
            RodElement().compute(self._params(**{"Number of Control Rods": 0}))

    def test_rod_shorter_than_absorber_rejected(self):
        with _DensityPatch(self.DENSITIES), self.assertRaises(ValueError):
            RodElement().compute(
                self._params(**{"Control Rod Length": 10.0,
                                "Control Rod Absorber Length": 50.0})
            )


class ResolveControlElementRodPathTest(unittest.TestCase):
    def test_sfr_routes_to_rod_and_sizes_it(self):
        params = {
            "reactor type": "SFR",
            "Number of Control Rods": 6,
            "Control Rod Radius": 2.0,
            "Control Rod Clad Thickness": 0.2,
            "Control Rod Absorber": "B4C_enriched",
            "Control Rod Clad": "SS304",
            "Active Height": 80.0,
            "Axial Reflector Thickness": 15.0,
        }
        with _DensityPatch({"B4C_enriched": 2.52, "SS304": 7.93}):
            strategy = resolve_control_element(params)
        self.assertIsInstance(strategy, RodElement)
        self.assertEqual(params["Control Element Type"], "rod")
        self.assertGreater(params["Control Element Mass"], 0.0)
        # Drum-specific keys must NOT be produced for a rod-controlled concept.
        self.assertNotIn("Control Drums Mass", params)


if __name__ == "__main__":
    unittest.main()
