# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
"""Control-element abstraction (control drums vs. control rods).

MOUSE's reactivity-control mass/geometry was previously synonymous with
*control drums*: ``core_design/drums.py`` assumes a rotating absorber-arc drum
embedded in the reflector, and the reflector/vessel geometry is even derived
from the drum placement. That model does not describe a sodium-cooled fast
reactor or an SRE-like concept, which use *control rods* inserted into the core
with drive mechanisms above the vessel.

This module introduces a small strategy interface, :class:`ControlElement`, with
two implementations:

* :class:`DrumElement` wraps the existing, validated drum calculation
  (``drums.calculate_drums_volumes_and_masses``) unchanged, so the LTMR/GCMR/
  HPMR numbers are byte-for-byte identical after the refactor.
* :class:`RodElement` is a new, self-contained control-rod model (absorber slug
  + clad + per-rod drive mechanism) that does not depend on any drum geometry.

Both implementations publish a common, element-agnostic vocabulary on the
params dict::

    Control Element Type   'drum' | 'rod'
    Control Element Count  number of elements
    Control Element Mass   total mass (kg), including drive mechanisms for rods

so downstream code can read reactivity-control mass without knowing whether the
concept uses drums or rods. ``DrumElement`` additionally keeps writing the
legacy ``Control Drums Mass`` / ``Drum Count`` / ``All Drums Area`` keys that
the cost layer consumes today.

The concrete strategy is chosen from the reactor capability registry, so adding
a rod-controlled concept needs no new ``if reactor type == ...`` branch here.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from core_design.openmc_materials_database import collect_materials_data
from core_design.drums import calculate_drums_volumes_and_masses
from core_design.reactor_registry import DRUM, ROD, get_capabilities, is_registered


class ControlElement(ABC):
    """Strategy interface for reactivity-control elements."""

    #: Registry value this strategy implements ('drum' or 'rod').
    kind: str = "control element"

    @abstractmethod
    def compute(self, params: dict) -> None:
        """Resolve geometry and populate mass keys on ``params`` in place."""

    @staticmethod
    def _publish_generic(params: dict, kind: str, count, mass: float) -> None:
        """Write the element-agnostic keys shared by every implementation."""
        params["Control Element Type"] = kind
        params["Control Element Count"] = count
        params["Control Element Mass"] = mass


class DrumElement(ControlElement):
    """Control drums embedded in the reflector (LTMR / GCMR / HPMR).

    Delegates entirely to the existing drum sizing/mass routine so behaviour is
    unchanged, then mirrors the result onto the generic control-element keys.
    """

    kind = DRUM

    def compute(self, params: dict) -> None:
        calculate_drums_volumes_and_masses(params)
        self._publish_generic(
            params,
            kind=DRUM,
            count=params.get("Drum Count"),
            mass=params["Control Drums Mass"],
        )


class RodElement(ControlElement):
    """Control rods inserted into the core, with drive mechanisms.

    Geometry per rod: a cylindrical absorber slug of radius
    ``r_abs = Control Rod Radius - Control Rod Clad Thickness`` over
    ``Control Rod Absorber Length``, inside a clad annulus over the full
    ``Control Rod Length``. A per-rod ``Control Rod Drive Mass`` accounts for
    the drive mechanism / extension above the absorber.

    Required params:
        Number of Control Rods, Control Rod Radius, Control Rod Clad Thickness,
        Control Rod Absorber, Control Rod Clad
    Optional params (sensible defaults derived from core geometry):
        Control Rod Absorber Length  (default: Active Height)
        Control Rod Length           (default: absorber length + 2*Axial
                                      Reflector Thickness when available)
        Control Rod Drive Mass       (kg per rod, default 0.0)

    Outputs (in addition to the generic keys):
        Control Rod Absorber Mass, Control Rod Clad Mass, Control Rods Mass,
        Control Rod Drive Mass Total, All Rods Volume, All Rods Area
    """

    kind = ROD

    def compute(self, params: dict) -> None:
        n_rods = int(params["Number of Control Rods"])
        if n_rods <= 0:
            raise ValueError(
                f"Number of Control Rods must be positive, got {n_rods}."
            )

        r_outer = float(params["Control Rod Radius"])
        clad_t = float(params["Control Rod Clad Thickness"])
        if clad_t <= 0:
            raise ValueError(
                f"Control Rod Clad Thickness must be positive, got {clad_t}."
            )
        r_abs = r_outer - clad_t
        if r_abs <= 0:
            raise ValueError(
                f"Control Rod Clad Thickness ({clad_t} cm) must be smaller than "
                f"Control Rod Radius ({r_outer} cm); absorber radius would be "
                f"{r_abs} cm."
            )

        absorber_length = float(
            params.get("Control Rod Absorber Length", params["Active Height"])
        )
        default_full_length = absorber_length + 2.0 * float(
            params.get("Axial Reflector Thickness", 0.0)
        )
        rod_length = float(params.get("Control Rod Length", default_full_length))
        if rod_length < absorber_length:
            raise ValueError(
                f"Control Rod Length ({rod_length} cm) cannot be shorter than "
                f"the absorber length ({absorber_length} cm)."
            )

        drive_mass_per_rod = float(params.get("Control Rod Drive Mass", 0.0))

        materials = collect_materials_data(params)
        absorber_density = materials[params["Control Rod Absorber"]].density  # g/cm^3
        clad_density = materials[params["Control Rod Clad"]].density          # g/cm^3

        absorber_vol = np.pi * r_abs * r_abs * absorber_length * n_rods
        clad_vol = np.pi * (r_outer * r_outer - r_abs * r_abs) * rod_length * n_rods

        # density in g/cm^3, volume in cm^3 -> /1000 gives kg (matches drums.py)
        absorber_mass = absorber_vol * absorber_density / 1000.0
        clad_mass = clad_vol * clad_density / 1000.0
        drive_mass_total = n_rods * drive_mass_per_rod

        rods_mass = absorber_mass + clad_mass
        total_mass = rods_mass + drive_mass_total

        params["Control Rod Absorber Mass"] = absorber_mass
        params["Control Rod Clad Mass"] = clad_mass
        params["Control Rods Mass"] = rods_mass
        params["Control Rod Drive Mass Total"] = drive_mass_total
        params["All Rods Volume"] = np.pi * r_outer * r_outer * rod_length * n_rods
        params["All Rods Area"] = np.pi * r_outer * r_outer * n_rods

        self._publish_generic(params, kind=ROD, count=n_rods, mass=total_mass)


# Strategy lookup keyed by the registry's control-element value.
_STRATEGIES = {
    DRUM: DrumElement,
    ROD: RodElement,
}


def _resolve_kind(params: dict) -> str:
    """Determine which control-element strategy to use.

    Precedence: an explicit ``params['Control Element Type']`` wins; otherwise
    the reactor capability registry decides; otherwise we fall back to drums for
    backward compatibility with callers that predate the registry.
    """
    explicit = params.get("Control Element Type")
    if explicit in _STRATEGIES:
        return explicit

    rtype = params.get("reactor type")
    if rtype and is_registered(rtype):
        return get_capabilities(rtype).control_element

    return DRUM


def resolve_control_element(params: dict) -> ControlElement:
    """Size the reactivity-control elements for ``params`` and return the
    strategy used.

    This is the capability-aware replacement for calling
    ``calculate_drums_volumes_and_masses`` directly: drum-controlled concepts
    behave exactly as before, while rod-controlled concepts route to
    :class:`RodElement`.
    """
    kind = _resolve_kind(params)
    strategy = _STRATEGIES[kind]()
    strategy.compute(params)
    return strategy
