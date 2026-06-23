"""Tests for the reactor capability registry."""

from __future__ import annotations

import unittest

from core_design.reactor_registry import (
    ReactorCapabilities,
    REACTOR_CAPABILITIES,
    get_capabilities,
    is_registered,
    implemented_reactor_types,
    THERMAL, FAST, GRAPHITE, METAL_HYDRIDE,
    GAS, LIQUID_METAL, HEAT_PIPE,
    DRUM, ROD, CONCENTRIC_RVACS, POOL, LOOP,
    TRISO, CYLINDRICAL_PIN, METALLIC,
)


class RegistryStructureTest(unittest.TestCase):
    def test_three_production_reactors_are_implemented(self):
        self.assertEqual(
            {"LTMR", "GCMR", "HPMR"}, set(implemented_reactor_types())
        )

    def test_sfr_and_sre_catalogued_but_not_implemented(self):
        for rt in ("SFR", "SRE"):
            self.assertTrue(is_registered(rt))
            self.assertFalse(get_capabilities(rt).implemented)

    def test_get_capabilities_unknown_raises_with_known_list(self):
        with self.assertRaises(KeyError) as ctx:
            get_capabilities("PWR")
        self.assertIn("LTMR", str(ctx.exception))

    def test_capabilities_are_frozen(self):
        caps = get_capabilities("LTMR")
        with self.assertRaises(Exception):
            caps.spectrum = FAST  # frozen dataclass


class CapabilityAxesTest(unittest.TestCase):
    def test_ltmr_axes(self):
        c = get_capabilities("LTMR")
        self.assertEqual(c.spectrum, THERMAL)
        self.assertEqual(c.moderator, METAL_HYDRIDE)
        self.assertEqual(c.coolant_class, LIQUID_METAL)
        self.assertEqual(c.control_element, DRUM)
        self.assertEqual(c.vessel_architecture, CONCENTRIC_RVACS)
        self.assertEqual(c.fuel_form, CYLINDRICAL_PIN)
        self.assertFalse(c.pressurized)

    def test_gcmr_is_pressurized_gas_triso(self):
        c = get_capabilities("GCMR")
        self.assertEqual(c.coolant_class, GAS)
        self.assertEqual(c.fuel_form, TRISO)
        self.assertTrue(c.pressurized)

    def test_hpmr_is_heat_pipe(self):
        self.assertEqual(get_capabilities("HPMR").coolant_class, HEAT_PIPE)

    def test_sfr_is_fast_rod_pool_metallic(self):
        c = get_capabilities("SFR")
        self.assertTrue(c.is_fast)
        self.assertFalse(c.is_moderated)
        self.assertTrue(c.uses_rods)
        self.assertFalse(c.uses_drums)
        self.assertEqual(c.vessel_architecture, POOL)
        self.assertEqual(c.fuel_form, METALLIC)
        self.assertEqual(c.coolant_class, LIQUID_METAL)

    def test_sre_is_thermal_graphite_loop_rod(self):
        c = get_capabilities("SRE")
        self.assertEqual(c.spectrum, THERMAL)
        self.assertEqual(c.moderator, GRAPHITE)
        self.assertTrue(c.is_moderated)
        self.assertTrue(c.uses_rods)
        self.assertEqual(c.vessel_architecture, LOOP)
        self.assertEqual(c.coolant_class, LIQUID_METAL)
        self.assertEqual(c.fuel_form, METALLIC)


class CapabilityValidationTest(unittest.TestCase):
    def test_invalid_axis_value_rejected(self):
        with self.assertRaises(ValueError):
            ReactorCapabilities(
                reactor_type="BAD",
                spectrum="luke-warm",  # invalid
                moderator=GRAPHITE,
                coolant_class=GAS,
                control_element=DRUM,
                vessel_architecture=CONCENTRIC_RVACS,
                fuel_form=TRISO,
            )

    def test_none_moderator_is_valid_for_fast_reactor(self):
        caps = ReactorCapabilities(
            reactor_type="X",
            spectrum=FAST,
            moderator=None,
            coolant_class=LIQUID_METAL,
            control_element=ROD,
            vessel_architecture=POOL,
            fuel_form=METALLIC,
        )
        self.assertIsNone(caps.moderator)
        self.assertFalse(caps.is_moderated)

    def test_every_registry_entry_validates(self):
        # Construction already validates; this asserts the registry is internally
        # consistent and uses only the documented axis vocabulary.
        for rt, caps in REACTOR_CAPABILITIES.items():
            self.assertEqual(rt, caps.reactor_type)


if __name__ == "__main__":
    unittest.main()
