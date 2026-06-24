# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

from core_design.reactor_registry import (
    is_registered,
    get_capabilities,
    CONCENTRIC_RVACS,
)
from reactor_engineering_evaluation.vessel_stack import (
    build_vessel_layers,
    compute_vessel_stack,
)


# Legacy per-layer output keys, keyed by the vessel-stack layer name. Kept so
# the cost layer and webapp continue to read the same params they always have.
_LEGACY_MASS_KEYS = {
    "inner": "Vessel Mass",
    "guard": "Guard Vessel Mass",
    "cooling": "Cooling Vessel Mass",
    "intake": "Intake Vessel Mass",
}
_LEGACY_RADIUS_KEYS = {
    "guard": "Guard Vessel Radius",
    "cooling": "Cooling Vessel Radius",
    "intake": "Intake Vessel Radius",
}


def _resolve_vessel_architecture(params):
    """Pick (architecture, pressurized) from the reactor capability registry.

    Falls back to the historical behaviour (concentric RVACS stack, with only
    the GCMR using pressurized full-ellipsoid heads) for any reactor-type string
    that is not in the registry, so legacy example scripts keep working.
    """
    rtype = params.get("reactor type")
    if rtype and is_registered(rtype):
        caps = get_capabilities(rtype)
        return caps.vessel_architecture, caps.pressurized
    return CONCENTRIC_RVACS, (rtype == "GCMR")


# Vessel Calcs
def vessels_specs(params):
    """Size the reactor vessel system and populate vessel mass/geometry params.

    The vessel system is described as an ordered stack of shells chosen by the
    reactor's ``vessel_architecture`` capability (concentric RVACS / pool /
    loop) and integrated by the generic vessel-stack walker. For the historical
    concentric architecture this reproduces the previous four-vessel result
    exactly, and the same legacy params (``Vessel Mass``, ``Guard Vessel
    Mass``, ..., ``Total Vessels Mass``, ``RPV Outer Radius`` for the GCMR) are
    still written.
    """
    # Shared cylindrical height for every shell in the stack (unchanged).
    vessel_height = (
        params["Active Height"]
        + 2 * params["Axial Reflector Thickness"]
        + params["Vessel Lower Plenum Height"]
        + params["Vessel Upper Plenum Height"]
        + params["Vessel Upper Gas Gap"]
    )

    architecture, pressurized = _resolve_vessel_architecture(params)
    layers = build_vessel_layers(params, architecture, pressurized)
    stack = compute_vessel_stack(layers, vessel_height)

    # --- Stack-level outputs (architecture-independent) ---
    params["Vessel Height"] = vessel_height
    params["Vessels Total Radius"] = stack.outer_radius
    params["Vessels Total Height"] = stack.total_height
    params["Total Vessels Mass"] = stack.total_mass
    # Generic breakdown for any number/kind of layers.
    params["Vessel Stack"] = [
        {
            "name": r.layer.name,
            "role": r.layer.role,
            "radius": r.layer.radius,
            "thickness": r.layer.thickness,
            "material": r.layer.material,
            "mass": r.mass,
        }
        for r in stack.layers
    ]

    # --- Legacy named outputs (written when the corresponding layer exists) ---
    present = {r.layer.name for r in stack.layers}
    for name, key in _LEGACY_MASS_KEYS.items():
        if name in present:
            params[key] = stack.mass_by_name(name)
    for name, key in _LEGACY_RADIUS_KEYS.items():
        if name in present:
            params[key] = stack.by_name(name).layer.radius

    # --- GCMR-style RPV reporting for a pressurized concentric stack ---
    if pressurized and architecture == CONCENTRIC_RVACS and "guard" in present:
        guard = stack.by_name("guard").layer
        params["RPV Outer Radius"] = guard.radius + guard.thickness
        params["RPV Outer Height"] = (
            vessel_height
            + 2 * params["Gap Between Vessel And Guard Vessel"]
            + 2 * params["Guard Vessel Thickness"]
            + 2 * params["Vessel Bottom Depth"]
            + 2 * params["Vessel Thickness"]
        )
