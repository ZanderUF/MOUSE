"""Tests for the SRE-like reactor build (sodium-cooled, graphite-moderated,
metallic-fuel, control-rod, loop-type concept).

Exercises the capability-driven abstractions end to end: the SRE registry entry
routes reactivity control to RodElement and the vessel system to the loop
architecture, and the builder produces a params dict that flows through the full
bottom-up cost estimate.
"""

from __future__ import annotations

import math
import unittest

import reactor_test_support as support

support.install()  # stub openmc/watts before importing the build/cost modules

from reactor_engineering_evaluation.fuel_calcs import (
    burnup_limited_fuel_lifetime_days,
)
from reactor_engineering_evaluation.tools import material_specific_heat


class BurnupLifetimeHelperTest(unittest.TestCase):
    def test_basic_formula(self):
        # 3 MWd/kgHM * 2000 kg / 20 MWt = 300 days
        self.assertAlmostEqual(
            burnup_limited_fuel_lifetime_days(20.0, 2000.0, 3.0), 300.0)

    def test_rejects_nonpositive_inputs(self):
        for bad in [
            (0.0, 2000.0, 3.0), (20.0, 0.0, 3.0), (20.0, 2000.0, 0.0),
        ]:
            with self.assertRaises(ValueError):
                burnup_limited_fuel_lifetime_days(*bad)


class SodiumPropertyTest(unittest.TestCase):
    def test_sodium_specific_heat_available(self):
        self.assertEqual(material_specific_heat("Na"), 1270.0)
        self.assertEqual(material_specific_heat("sodium"), 1270.0)


class SreBuildTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = support.build_case("SRE")

    def test_identity_and_capabilities(self):
        p = self.params
        self.assertEqual(p["reactor type"], "SRE")
        self.assertEqual(p["Coolant"], "Na")
        self.assertEqual(p["Fuel"], "U_met")
        self.assertEqual(p["Moderator"], "Graphite")

    def test_uses_control_rods_not_drums(self):
        p = self.params
        self.assertEqual(p["Control Element Type"], "rod")
        self.assertEqual(p["Control Element Count"], p["Number of Control Rods"])
        self.assertGreater(p["Control Rods Mass"], 0)
        self.assertIn("Control Rod Absorber Mass", p)
        # Rods have no reflector face, so no drum-reflector mass is produced.
        self.assertNotIn("Control Drum Reflector Mass", p)

    def test_uses_loop_vessel_architecture(self):
        stack = self.params["Vessel Stack"]
        self.assertEqual([layer["name"] for layer in stack], ["inner"])
        self.assertAlmostEqual(self.params["Total Vessels Mass"],
                               self.params["Vessel Mass"])

    def test_metallic_fuel_heavy_metal_inventory(self):
        # Historical SRE Core I held roughly 2.7-2.9 tonnes of uranium metal.
        self.assertTrue(2500 < self.params["Uranium Mass"] < 3200,
                        self.params["Uranium Mass"])

    def test_burnup_limited_fuel_lifetime(self):
        p = self.params
        expected = burnup_limited_fuel_lifetime_days(
            p["Power MWt"], p["Uranium Mass"],
            p["Target Discharge Burnup MWd/kgHM"])
        self.assertAlmostEqual(p["Fuel Lifetime"], round(expected), delta=1)
        self.assertGreaterEqual(p["Fuel Lifetime"], 90)

    def test_graphite_masses_positive(self):
        for key in ("Moderator Mass", "Radial Reflector Mass", "Axial Reflector Mass"):
            self.assertGreater(self.params[key], 0, key)

    def test_loop_heat_transport_sized(self):
        p = self.params
        self.assertEqual(p["Primary Pump"], "Yes")
        self.assertEqual(p["Secondary Pump"], "Yes")  # intermediate sodium loop
        self.assertGreater(p["Primary HX Mass"], 0)
        self.assertGreater(p["Primary Pump Mechanical Power"], 0)
        self.assertGreater(p["Secondary Pump Mechanical Power"], 0)


class SreCostEstimateTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        support.install()
        from webapp.estimate_service import EstimateInputs, run_estimate

        base = dict(
            interest_rate=0.05, discount_rate=0.05, construction_duration=12,
            debt_to_equity=0.5, operation_mode="Remotely Monitored",
            emergency_shutdowns=2.0, startup_duration=21,
            startup_duration_refueling=30, tax_credit_type="ITC",
            tax_credit_value=0.30, plant_lifetime=60, tax_credit_units=10,
        )
        cls.result = run_estimate(EstimateInputs(
            reactor_type="SRE", power_mwt=20, enrichment=0.0278,
            active_height=180, **base))

    def test_cost_tables_non_empty(self):
        self.assertFalse(self.result.display_df.empty)
        self.assertFalse(self.result.enriched_df.empty)
        self.assertEqual(self.result.params["reactor type"], "SRE")

    def test_lcoe_finite_and_positive(self):
        ed = self.result.enriched_df
        for col in ("FOAK LCOE", "NOAK LCOE"):
            self.assertIn(col, ed.columns)
            total = float(ed[col].sum())
            self.assertTrue(math.isfinite(total))
            self.assertGreater(total, 0.0)
            self.assertTrue((ed[col].dropna() >= 0).all())

    def test_capacity_factor_in_range(self):
        cf = self.result.params["Capacity Factor"]
        self.assertTrue(0 < cf <= 1, cf)


if __name__ == "__main__":
    unittest.main()
