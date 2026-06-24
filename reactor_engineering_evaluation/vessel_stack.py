# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""Vessel-stack model.

The original ``vessels_specs`` hard-coded exactly four concentric vessels
(inner / guard / cooling / intake) with the GCMR special-cased for
full-ellipsoid heads. A pool-type sodium fast reactor (one large primary vessel
plus a guard vessel, no RVACS air shells) or a loop-type SRE-like concept does
not fit that fixed four-shell template.

This module decomposes the vessel system into an ordered list of
:class:`VesselLayer` records and a single generic walker,
:func:`compute_vessel_stack`, that sizes any number of layers using the same
ellipsoidal-head + cylindrical-shell formulas as before. Each *architecture* is
just a function that returns the appropriate layer list:

* :func:`concentric_rvacs_layers` reproduces the historical four-vessel stack
  bit-for-bit (including its radius/bottom-depth chaining quirk and the
  pressurized full-ellipsoid heads).
* :func:`pool_layers` and :func:`loop_layers` describe single-primary-vessel
  architectures for the sodium concepts.

The walker is deliberately free of chaining assumptions: the architecture
builder computes each layer's explicit inner radius, wall thickness, head depth,
head factor, and material, and the walker only integrates them. New
architectures therefore need no change to the walker.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Tuple

from reactor_engineering_evaluation.tools import (
    ellipsoid_shell,
    circle_area,
    materials_densities,
)


# Head geometry factors: a closed pressure vessel models both the bottom head
# and (effectively) the top closure, while an open-top pool/plenum models only
# the bottom head.
HEAD_FULL = 1.0   # full ellipsoidal head (closed pressure vessel)
HEAD_HALF = 0.5   # half ellipsoidal head (open-top pool / plenum)
HEAD_NONE = 0.0   # no ellipsoidal head contribution


@dataclass
class VesselLayer:
    """A single cylindrical shell with an ellipsoidal bottom head.

    All dimensions in cm. ``radius`` is the *inner* radius of the shell.
    """

    name: str
    radius: float
    thickness: float
    bottom_depth: float
    material: str
    head_factor: float = HEAD_HALF
    role: str = "structural"  # 'structural' | 'cooling' | 'intake' | ...


@dataclass
class VesselLayerResult:
    layer: VesselLayer
    volume: float  # cm^3
    mass: float    # kg


@dataclass
class VesselStackResult:
    layers: List[VesselLayerResult]
    total_mass: float        # kg
    outer_radius: float      # cm (outermost shell outer radius)
    total_height: float      # cm (outer head depth + cylindrical height)
    cylindrical_height: float  # cm (shared shell height used for every layer)

    def mass_by_name(self, name: str) -> float:
        for r in self.layers:
            if r.layer.name == name:
                return r.mass
        raise KeyError(f"No vessel layer named {name!r} in stack.")

    def by_name(self, name: str) -> VesselLayerResult:
        for r in self.layers:
            if r.layer.name == name:
                return r
        raise KeyError(f"No vessel layer named {name!r} in stack.")


def _layer_volume(layer: VesselLayer, cylindrical_height: float) -> float:
    """Ellipsoidal-head shell volume + cylindrical-wall shell volume (cm^3).

    Identical to the per-vessel expressions in the original ``vessels_specs``:
    ``ellipsoid_shell(r, r, depth) * head_factor * thickness`` for the head and
    ``(circle_area(r + t) - circle_area(r)) * cylindrical_height`` for the wall.
    """
    head = (
        ellipsoid_shell(layer.radius, layer.radius, layer.bottom_depth)
        * layer.head_factor
        * layer.thickness
    )
    wall = (
        circle_area(layer.radius + layer.thickness) - circle_area(layer.radius)
    ) * cylindrical_height
    return head + wall


def compute_vessel_stack(
    layers: List[VesselLayer],
    cylindrical_height: float,
    density_fn: Callable[[str], float] = materials_densities,
) -> VesselStackResult:
    """Integrate masses for an ordered list of vessel layers.

    ``density_fn`` maps a material name to a density in g/cm^3 (default: the
    engineering-layer ``materials_densities`` table). Masses are returned in kg
    (volume[cm^3] * density[g/cm^3] / 1000), matching the existing convention.

    The outermost layer (last in the list) defines the stack's overall outer
    radius and height.
    """
    if not layers:
        raise ValueError("compute_vessel_stack requires at least one layer.")

    results: List[VesselLayerResult] = []
    total_mass = 0.0
    for layer in layers:
        vol = _layer_volume(layer, cylindrical_height)
        mass = vol * density_fn(layer.material) / 1000.0
        results.append(VesselLayerResult(layer=layer, volume=vol, mass=mass))
        total_mass += mass

    outer = layers[-1]
    outer_radius = outer.radius + outer.thickness
    total_height = outer.bottom_depth + cylindrical_height

    return VesselStackResult(
        layers=results,
        total_mass=total_mass,
        outer_radius=outer_radius,
        total_height=total_height,
        cylindrical_height=cylindrical_height,
    )


# --------------------------------------------------------------------------- #
# Architecture layer builders
# --------------------------------------------------------------------------- #

def concentric_rvacs_layers(params: dict, pressurized: bool) -> List[VesselLayer]:
    """Reproduce the historical inner/guard/cooling/intake concentric stack.

    The radius and bottom-depth chaining intentionally matches the original
    ``vessels_specs`` exactly, including the asymmetry where the cooling and
    intake *radii* do not add the preceding wall thickness while the *bottom
    depths* do. When ``pressurized`` the inner and guard vessels use full
    ellipsoidal heads (the GCMR RPV); otherwise half heads. Cooling and intake
    shells always use half heads.
    """
    R = params["Vessel Radius"]
    t = params["Vessel Thickness"]
    bd = params["Vessel Bottom Depth"]

    gap_vg = params["Gap Between Vessel And Guard Vessel"]
    Gt = params["Guard Vessel Thickness"]
    gap_gc = params["Gap Between Guard Vessel And Cooling Vessel"]
    Ct = params["Cooling Vessel Thickness"]
    gap_ci = params["Gap Between Cooling Vessel And Intake Vessel"]
    It = params["Intake Vessel Thickness"]

    structural_head = HEAD_FULL if pressurized else HEAD_HALF

    guard_radius = R + t + gap_vg
    guard_bottom = bd + t + gap_vg
    # NOTE: original does NOT add the preceding wall thickness to the radius
    # here, but DOES add it to the bottom depth. Preserved deliberately.
    cooling_radius = guard_radius + gap_gc
    cooling_bottom = guard_bottom + Gt + gap_gc
    intake_radius = cooling_radius + gap_ci
    intake_bottom = cooling_bottom + Ct + gap_ci

    return [
        VesselLayer("inner", R, t, bd, params["Vessel Material"],
                    structural_head, "structural"),
        VesselLayer("guard", guard_radius, Gt, guard_bottom,
                    params["Guard Vessel Material"], structural_head, "structural"),
        VesselLayer("cooling", cooling_radius, Ct, cooling_bottom,
                    params["Cooling Vessel Material"], HEAD_HALF, "cooling"),
        VesselLayer("intake", intake_radius, It, intake_bottom,
                    params["Intake Vessel Material"], HEAD_HALF, "intake"),
    ]


def pool_layers(params: dict) -> List[VesselLayer]:
    """Pool-type architecture: one primary vessel plus a guard vessel.

    Suitable for a sodium-cooled fast microreactor where the core, primary
    pumps, and intermediate heat exchanger sit inside a single large
    low-pressure vessel, with a guard vessel to contain a primary-coolant leak.
    No RVACS air shells. Uses half ellipsoidal heads (open-top pool with a cover
    gas above the free sodium surface).
    """
    R = params["Vessel Radius"]
    t = params["Vessel Thickness"]
    bd = params["Vessel Bottom Depth"]
    gap_vg = params["Gap Between Vessel And Guard Vessel"]
    Gt = params["Guard Vessel Thickness"]

    guard_radius = R + t + gap_vg
    guard_bottom = bd + t + gap_vg

    return [
        VesselLayer("inner", R, t, bd, params["Vessel Material"],
                    HEAD_HALF, "structural"),
        VesselLayer("guard", guard_radius, Gt, guard_bottom,
                    params["Guard Vessel Material"], HEAD_HALF, "structural"),
    ]


def loop_layers(params: dict) -> List[VesselLayer]:
    """Loop-type architecture: a single primary vessel (no guard, no RVACS).

    Suitable for an SRE-like concept where the primary coolant circulates
    through external loops and heat exchangers rather than a pool. Only the
    primary vessel mass is modelled here; external loop components are accounted
    for elsewhere in the balance-of-plant.
    """
    R = params["Vessel Radius"]
    t = params["Vessel Thickness"]
    bd = params["Vessel Bottom Depth"]

    return [
        VesselLayer("inner", R, t, bd, params["Vessel Material"],
                    HEAD_HALF, "structural"),
    ]


# Architecture name -> layer builder. ``concentric_rvacs`` takes the
# ``pressurized`` flag, handled in build_vessel_layers below.
def build_vessel_layers(params: dict, architecture: str,
                        pressurized: bool) -> List[VesselLayer]:
    """Return the vessel layer list for the named architecture."""
    if architecture == "concentric_rvacs":
        return concentric_rvacs_layers(params, pressurized)
    if architecture == "pool":
        return pool_layers(params)
    if architecture == "loop":
        return loop_layers(params)
    raise ValueError(
        f"Unknown vessel architecture {architecture!r}; expected one of "
        "'concentric_rvacs', 'pool', 'loop'."
    )
