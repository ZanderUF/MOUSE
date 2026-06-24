"""Tests for the generic vessel-stack model and the refactored vessels_specs."""

from __future__ import annotations

import unittest

from reactor_engineering_evaluation.tools import (
    ellipsoid_shell,
    circle_area,
    materials_densities,
)
from reactor_engineering_evaluation.vessel_stack import (
    VesselLayer,
    compute_vessel_stack,
    concentric_rvacs_layers,
    pool_layers,
    loop_layers,
    build_vessel_layers,
    HEAD_FULL,
    HEAD_HALF,
)
from reactor_engineering_evaluation.vessels_calcs import vessels_specs


def _vessel_inputs(**overrides):
    """A complete set of vessel input params (concentric RVACS template)."""
    p = {
        "Active Height": 100.0,
        "Axial Reflector Thickness": 20.0,
        "Vessel Lower Plenum Height": 50.0,
        "Vessel Upper Plenum Height": 47.0,
        "Vessel Upper Gas Gap": 0.0,
        "Vessel Radius": 60.0,
        "Vessel Thickness": 2.0,
        "Vessel Bottom Depth": 32.0,
        "Vessel Material": "stainless_steel",
        "Gap Between Vessel And Guard Vessel": 5.0,
        "Guard Vessel Thickness": 1.0,
        "Guard Vessel Material": "stainless_steel",
        "Gap Between Guard Vessel And Cooling Vessel": 5.0,
        "Cooling Vessel Thickness": 0.5,
        "Cooling Vessel Material": "stainless_steel",
        "Gap Between Cooling Vessel And Intake Vessel": 5.0,
        "Intake Vessel Thickness": 0.5,
        "Intake Vessel Material": "stainless_steel",
    }
    p.update(overrides)
    return p


def _reference_concentric(params, pressurized):
    """Independent reimplementation of the ORIGINAL four-vessel formula.

    Mirrors the pre-refactor vessels_specs exactly so we can assert the new
    stack-based implementation produces identical masses.
    """
    vh = (params["Active Height"] + 2 * params["Axial Reflector Thickness"]
          + params["Vessel Lower Plenum Height"]
          + params["Vessel Upper Plenum Height"]
          + params["Vessel Upper Gas Gap"])
    R, t, bd = params["Vessel Radius"], params["Vessel Thickness"], params["Vessel Bottom Depth"]
    dens = materials_densities

    def wall(r, th):
        return (circle_area(r + th) - circle_area(r)) * vh

    if pressurized:
        vv = ellipsoid_shell(R, R, bd) * t + wall(R, t)
    else:
        vv = ellipsoid_shell(R, R, bd) / 2 * t + wall(R, t)
    vm = vv * dens(params["Vessel Material"]) / 1000

    gr = R + t + params["Gap Between Vessel And Guard Vessel"]
    gbd = bd + t + params["Gap Between Vessel And Guard Vessel"]
    Gt = params["Guard Vessel Thickness"]
    if pressurized:
        gv = ellipsoid_shell(gr, gr, gbd) * Gt + wall(gr, Gt)
    else:
        gv = ellipsoid_shell(gr, gr, gbd) / 2 * Gt + wall(gr, Gt)
    gm = gv * dens(params["Guard Vessel Material"]) / 1000

    cr = gr + params["Gap Between Guard Vessel And Cooling Vessel"]
    cbd = gbd + Gt + params["Gap Between Guard Vessel And Cooling Vessel"]
    Ct = params["Cooling Vessel Thickness"]
    cv = ellipsoid_shell(cr, cr, cbd) / 2 * Ct + wall(cr, Ct)
    cm = cv * dens(params["Cooling Vessel Material"]) / 1000

    ir = cr + params["Gap Between Cooling Vessel And Intake Vessel"]
    ibd = cbd + Ct + params["Gap Between Cooling Vessel And Intake Vessel"]
    It = params["Intake Vessel Thickness"]
    iv = ellipsoid_shell(ir, ir, ibd) / 2 * It + wall(ir, It)
    im = iv * dens(params["Intake Vessel Material"]) / 1000

    return {
        "Vessel Mass": vm, "Guard Vessel Mass": gm,
        "Cooling Vessel Mass": cm, "Intake Vessel Mass": im,
        "Total Vessels Mass": vm + gm + cm + im,
        "Guard Vessel Radius": gr, "Cooling Vessel Radius": cr,
        "Intake Vessel Radius": ir,
        "Vessels Total Radius": ir + It, "Vessels Total Height": ibd + vh,
        "Vessel Height": vh,
    }


class WalkerTest(unittest.TestCase):
    def test_single_layer_volume_and_mass(self):
        layer = VesselLayer("inner", radius=100.0, thickness=2.0,
                            bottom_depth=30.0, material="steel",
                            head_factor=HEAD_HALF)
        height = 200.0
        result = compute_vessel_stack([layer], height, density_fn=lambda m: 8.0)

        exp_head = ellipsoid_shell(100, 100, 30) * 0.5 * 2.0
        exp_wall = (circle_area(102) - circle_area(100)) * height
        exp_vol = exp_head + exp_wall
        self.assertAlmostEqual(result.layers[0].volume, exp_vol)
        self.assertAlmostEqual(result.layers[0].mass, exp_vol * 8.0 / 1000)
        self.assertAlmostEqual(result.outer_radius, 102.0)
        self.assertAlmostEqual(result.total_height, 230.0)

    def test_empty_stack_rejected(self):
        with self.assertRaises(ValueError):
            compute_vessel_stack([], 100.0)


class LayerBuilderTest(unittest.TestCase):
    def test_concentric_has_four_named_layers(self):
        layers = concentric_rvacs_layers(_vessel_inputs(), pressurized=False)
        self.assertEqual([l.name for l in layers],
                         ["inner", "guard", "cooling", "intake"])

    def test_concentric_chaining_matches_original(self):
        p = _vessel_inputs()
        layers = {l.name: l for l in concentric_rvacs_layers(p, pressurized=False)}
        gr = p["Vessel Radius"] + p["Vessel Thickness"] + p["Gap Between Vessel And Guard Vessel"]
        self.assertAlmostEqual(layers["guard"].radius, gr)
        # cooling/intake radii deliberately do NOT add preceding wall thickness
        self.assertAlmostEqual(layers["cooling"].radius,
                               gr + p["Gap Between Guard Vessel And Cooling Vessel"])

    def test_pressurized_uses_full_heads_for_structural_shells(self):
        layers = {l.name: l for l in concentric_rvacs_layers(_vessel_inputs(), pressurized=True)}
        self.assertEqual(layers["inner"].head_factor, HEAD_FULL)
        self.assertEqual(layers["guard"].head_factor, HEAD_FULL)
        self.assertEqual(layers["cooling"].head_factor, HEAD_HALF)

    def test_pool_has_inner_and_guard_only(self):
        layers = pool_layers(_vessel_inputs())
        self.assertEqual([l.name for l in layers], ["inner", "guard"])
        self.assertTrue(all(l.head_factor == HEAD_HALF for l in layers))

    def test_loop_has_single_vessel(self):
        layers = loop_layers(_vessel_inputs())
        self.assertEqual([l.name for l in layers], ["inner"])

    def test_build_vessel_layers_unknown_architecture(self):
        with self.assertRaises(ValueError):
            build_vessel_layers(_vessel_inputs(), "toroidal", False)


class VesselsSpecsConcentricRegressionTest(unittest.TestCase):
    """New vessels_specs must match the original four-vessel math exactly."""

    def _check(self, reactor_type, pressurized):
        params = _vessel_inputs(**{"reactor type": reactor_type})
        vessels_specs(params)
        ref = _reference_concentric(params, pressurized)
        for key, expected in ref.items():
            self.assertAlmostEqual(params[key], expected, places=6,
                                   msg=f"{reactor_type}: {key}")

    def test_non_pressurized_concentric_ltmr(self):
        self._check("LTMR", pressurized=False)

    def test_pressurized_concentric_gcmr_writes_rpv(self):
        self._check("GCMR", pressurized=True)
        params = _vessel_inputs(**{"reactor type": "GCMR"})
        vessels_specs(params)
        self.assertIn("RPV Outer Radius", params)
        self.assertIn("RPV Outer Height", params)


class VesselsSpecsAlternateArchitectureTest(unittest.TestCase):
    """The refactor enables pool/loop architectures via the same entry point."""

    def test_pool_architecture_sfr(self):
        # SFR is registered with a pool vessel architecture.
        params = _vessel_inputs(**{"reactor type": "SFR"})
        vessels_specs(params)

        self.assertEqual(len(params["Vessel Stack"]), 2)
        self.assertIn("Vessel Mass", params)
        self.assertIn("Guard Vessel Mass", params)
        # No RVACS shells and no pressurized RPV for a pool reactor.
        self.assertNotIn("Cooling Vessel Mass", params)
        self.assertNotIn("Intake Vessel Mass", params)
        self.assertNotIn("RPV Outer Radius", params)
        self.assertAlmostEqual(
            params["Total Vessels Mass"],
            params["Vessel Mass"] + params["Guard Vessel Mass"],
        )

    def test_loop_architecture_sre(self):
        # SRE is registered with a loop vessel architecture (single vessel).
        params = _vessel_inputs(**{"reactor type": "SRE"})
        vessels_specs(params)

        self.assertEqual(len(params["Vessel Stack"]), 1)
        self.assertAlmostEqual(params["Total Vessels Mass"], params["Vessel Mass"])
        self.assertNotIn("Guard Vessel Mass", params)
        self.assertNotIn("Cooling Vessel Mass", params)


if __name__ == "__main__":
    unittest.main()
