# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""Reactor capability registry.

Historically MOUSE encoded the reactor concept as a bare string
('LTMR' / 'GCMR' / 'HPMR') that was dispatched through ``if/elif`` chains and
per-type builder functions scattered across the webapp, engineering, and cost
layers. Physics architecture was therefore bundled under a *name* rather than
described as a set of *capabilities*, which made it expensive to add a new
concept (e.g. a sodium-cooled fast microreactor that uses control rods instead
of drums and a pool-type vessel instead of the concentric RVACS stack).

This module makes the reactor's physical character a first-class, data-driven
description. Each concept is a :class:`ReactorCapabilities` record whose
independent axes (spectrum, moderator, coolant class, control-element type,
vessel architecture, fuel form) can be combined freely. Subsystems
(control-element sizing, vessel-stack construction, ...) read these
capabilities instead of branching on the reactor name, so a new concept plugs
in by adding one registry entry plus the parameter builder for it.

The registry is intentionally free of heavy imports (numpy/openmc/pandas) so it
can be imported from any layer without pulling in the simulation stack.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple


# --------------------------------------------------------------------------- #
# Capability vocabulary. Using module-level constants (rather than bare string
# literals at every call site) keeps the allowed values discoverable and
# typo-proof.
# --------------------------------------------------------------------------- #

# Neutron spectrum
THERMAL = "thermal"
FAST = "fast"
SPECTRA = (THERMAL, FAST)

# Moderator (None for a fast reactor with no dedicated moderator)
GRAPHITE = "graphite"
METAL_HYDRIDE = "metal_hydride"  # ZrH, YHx, ...
MODERATORS = (None, GRAPHITE, METAL_HYDRIDE)

# Coolant class drives the balance-of-plant / power-cycle strategy
GAS = "gas"
LIQUID_METAL = "liquid_metal"
HEAT_PIPE = "heat_pipe"
COOLANT_CLASSES = (GAS, LIQUID_METAL, HEAT_PIPE)

# Reactivity-control element
DRUM = "drum"
ROD = "rod"
CONTROL_ELEMENTS = (DRUM, ROD)

# Vessel architecture
CONCENTRIC_RVACS = "concentric_rvacs"  # nested inner/guard/cooling/intake shells
POOL = "pool"                          # single large primary vessel + guard
LOOP = "loop"                          # primary vessel + external loop vessels
VESSEL_ARCHITECTURES = (CONCENTRIC_RVACS, POOL, LOOP)

# Fuel form
TRISO = "triso"
CYLINDRICAL_PIN = "cylindrical_pin"
METALLIC = "metallic"
FUEL_FORMS = (TRISO, CYLINDRICAL_PIN, METALLIC)

# Core lattice
HEX = "hex"
SQUARE = "square"
LATTICES = (HEX, SQUARE)


@dataclass(frozen=True)
class ReactorCapabilities:
    """Capability descriptor for a single reactor concept.

    Each field is an independent axis; subsystems read the axis they care about
    rather than the reactor name.
    """

    reactor_type: str
    spectrum: str
    moderator: Optional[str]
    coolant_class: str
    control_element: str
    vessel_architecture: str
    fuel_form: str
    lattice: str = HEX
    # Whether the primary vessel is a closed pressure boundary. Drives the
    # vessel head geometry (full vs half ellipsoid) and RPV reporting.
    pressurized: bool = False
    # True once a webapp/params builder exists for this concept. Catalogued but
    # unbuilt concepts (e.g. SFR, SRE) describe the target so the abstractions
    # can be exercised and the wiring point is obvious.
    implemented: bool = True
    description: str = ""

    def __post_init__(self):
        _validate("spectrum", self.spectrum, SPECTRA)
        _validate("moderator", self.moderator, MODERATORS)
        _validate("coolant_class", self.coolant_class, COOLANT_CLASSES)
        _validate("control_element", self.control_element, CONTROL_ELEMENTS)
        _validate("vessel_architecture", self.vessel_architecture, VESSEL_ARCHITECTURES)
        _validate("fuel_form", self.fuel_form, FUEL_FORMS)
        _validate("lattice", self.lattice, LATTICES)

    # Convenience predicates (read better at call sites than string compares)
    @property
    def is_fast(self) -> bool:
        return self.spectrum == FAST

    @property
    def is_moderated(self) -> bool:
        return self.moderator is not None

    @property
    def uses_drums(self) -> bool:
        return self.control_element == DRUM

    @property
    def uses_rods(self) -> bool:
        return self.control_element == ROD


def _validate(field_name: str, value, allowed: Tuple):
    if value not in allowed:
        raise ValueError(
            f"Invalid {field_name}={value!r}; expected one of {allowed}."
        )


# --------------------------------------------------------------------------- #
# Registry of reactor concepts.
#
# The three production concepts (LTMR/GCMR/HPMR) are flagged ``implemented``.
# SFR and SRE are catalogued targets (``implemented=False``): they document the
# capability combination the abstractions are being generalised toward and let
# tests exercise the rod/pool/loop paths, but ``build_params`` does not yet have
# a parameter builder for them.
# --------------------------------------------------------------------------- #

REACTOR_CAPABILITIES: Dict[str, ReactorCapabilities] = {
    "LTMR": ReactorCapabilities(
        reactor_type="LTMR",
        spectrum=THERMAL,
        moderator=METAL_HYDRIDE,
        coolant_class=LIQUID_METAL,
        control_element=DRUM,
        vessel_architecture=CONCENTRIC_RVACS,
        fuel_form=CYLINDRICAL_PIN,
        pressurized=False,
        implemented=True,
        description=(
            "Liquid-metal (NaK) thermal microreactor: ZrH-moderated, "
            "cylindrical pin fuel, control drums in a graphite reflector, "
            "concentric guard/cooling/intake (RVACS) vessels."
        ),
    ),
    "GCMR": ReactorCapabilities(
        reactor_type="GCMR",
        spectrum=THERMAL,
        moderator=GRAPHITE,
        coolant_class=GAS,
        control_element=DRUM,
        vessel_architecture=CONCENTRIC_RVACS,
        fuel_form=TRISO,
        pressurized=True,
        implemented=True,
        description=(
            "Gas-cooled (He) thermal microreactor: graphite-moderated, TRISO "
            "fuel, control drums, pressurized RPV with full-ellipsoid heads."
        ),
    ),
    "HPMR": ReactorCapabilities(
        reactor_type="HPMR",
        spectrum=THERMAL,
        moderator=GRAPHITE,
        coolant_class=HEAT_PIPE,
        control_element=DRUM,
        vessel_architecture=CONCENTRIC_RVACS,
        fuel_form=TRISO,
        pressurized=False,
        implemented=True,
        description=(
            "Heat-pipe thermal microreactor: graphite-moderated, TRISO fuel, "
            "control drums, passive (no primary pump) heat transport."
        ),
    ),
    # ------------------------------------------------------------------ #
    # Catalogued generalisation targets (no params builder yet).
    # ------------------------------------------------------------------ #
    "SFR": ReactorCapabilities(
        reactor_type="SFR",
        spectrum=FAST,
        moderator=None,
        coolant_class=LIQUID_METAL,
        control_element=ROD,
        vessel_architecture=POOL,
        fuel_form=METALLIC,
        pressurized=False,
        implemented=False,
        description=(
            "Sodium-cooled fast microreactor: unmoderated, metallic fuel, "
            "control rods with drive mechanisms, pool-type primary vessel."
        ),
    ),
    "SRE": ReactorCapabilities(
        reactor_type="SRE",
        spectrum=THERMAL,
        moderator=GRAPHITE,
        coolant_class=LIQUID_METAL,
        control_element=ROD,
        vessel_architecture=LOOP,
        fuel_form=METALLIC,
        pressurized=False,
        implemented=True,
        description=(
            "Sodium Reactor Experiment-like concept: graphite-moderated "
            "thermal spectrum, sodium coolant, metallic uranium fuel, control "
            "rods, loop-type primary system."
        ),
    ),
}


def get_capabilities(reactor_type: str) -> ReactorCapabilities:
    """Return the :class:`ReactorCapabilities` for ``reactor_type``.

    Raises ``KeyError`` with the list of known concepts if unknown.
    """
    try:
        return REACTOR_CAPABILITIES[reactor_type]
    except KeyError:
        known = ", ".join(sorted(REACTOR_CAPABILITIES))
        raise KeyError(
            f"Unknown reactor type {reactor_type!r}. Known concepts: {known}."
        ) from None


def is_registered(reactor_type: str) -> bool:
    return reactor_type in REACTOR_CAPABILITIES


def implemented_reactor_types() -> Tuple[str, ...]:
    """Reactor types that have a parameter builder (usable by build_params)."""
    return tuple(
        rt for rt, caps in REACTOR_CAPABILITIES.items() if caps.implemented
    )
