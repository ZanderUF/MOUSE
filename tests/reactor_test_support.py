"""Shared lightweight test harness for the reactor-config calculation path.

Mirrors the dependency strategy of tests/test_estimate_service.py: OpenMC /
WATTS / matplotlib are stubbed and material densities are sourced from
webapp/materials_densities.json, so the params-building and engineering
calculations can run without the full simulation stack installed.

Import this module *before* importing anything from webapp/core_design so the
runtime stubs are installed first.
"""

from __future__ import annotations

import json
import os
import sys
from unittest.mock import MagicMock


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEBAPP_DIR = os.path.join(REPO_ROOT, "webapp")
for _p in (REPO_ROOT, WEBAPP_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)


class _MaterialStub:
    def __init__(self, name=None, temperature=None):
        self.name = name
        self.temperature = temperature
        self.density = 0.0

    def set_density(self, units, value):
        self.density = value

    def add_nuclide(self, *args, **kwargs):
        pass

    def add_element(self, *args, **kwargs):
        pass

    def add_s_alpha_beta(self, *args, **kwargs):
        pass

    @staticmethod
    def mix_materials(materials, fractions, method, name=None):
        result = _MaterialStub(name=name)
        try:
            result.density = sum(m.density * f for m, f in zip(materials, fractions))
        except Exception:
            result.density = 0.0
        return result


class _MaterialsStub:
    def append(self, mat):
        pass

    def extend(self, mats):
        pass


class ThinMaterial:
    __slots__ = ("name", "density")

    def __init__(self, name, density):
        self.name = name
        self.density = density


_INSTALLED = False
_DENSITIES_RAW = None


def install():
    """Install runtime stubs and the JSON-backed material-density lookup.

    Idempotent: safe to call from multiple test modules.
    """
    global _INSTALLED, _DENSITIES_RAW
    if _INSTALLED:
        return _DENSITIES_RAW

    openmc_stub = MagicMock()
    openmc_stub.Material = _MaterialStub
    openmc_stub.Materials = _MaterialsStub
    for mod in ["openmc", "openmc.deplete", "openmc.mgxs"]:
        sys.modules[mod] = openmc_stub
    sys.modules["watts"] = MagicMock()

    mpl_stub = MagicMock()
    for mod in ["matplotlib", "matplotlib.pyplot", "matplotlib.patches",
                "matplotlib.colors"]:
        sys.modules[mod] = mpl_stub

    import core_design.openmc_materials_database as materials_db_mod

    with open(os.path.join(WEBAPP_DIR, "materials_densities.json")) as f:
        raw = json.load(f)

    thin_by_reactor = {
        rtype: {name: ThinMaterial(name, float(density))
                for name, density in mats.items()}
        for rtype, mats in raw.items()
    }

    def collect_materials_data(params):
        rtype = params.get("reactor type", "LTMR")
        return thin_by_reactor.get(rtype, thin_by_reactor["LTMR"])

    materials_db_mod.collect_materials_data = collect_materials_data

    _INSTALLED = True
    _DENSITIES_RAW = raw
    return raw


# Canonical design points (match tests/test_estimate_service.py).
REACTOR_BUILD_CASES = {
    "LTMR": dict(power_mwt=20, enrichment=0.1975,
                 n_rings_per_assembly=12, active_height=95),
    "GCMR": dict(power_mwt=15, enrichment=0.1975,
                 n_assembly_rings=6, n_core_rings=5, active_height=192),
    "HPMR": dict(power_mwt=5, enrichment=0.1975,
                 n_assembly_rings=6, n_core_rings=5, active_height=109),
}


def build_case(reactor_type):
    """Build a fully-populated params dict for a canonical design point."""
    install()
    from webapp.reactor_config import build_params

    case = REACTOR_BUILD_CASES[reactor_type]
    return build_params(
        reactor_type,
        case["power_mwt"],
        case["enrichment"],
        {},
        n_rings_per_assembly=case.get("n_rings_per_assembly"),
        active_height=case.get("active_height"),
        n_assembly_rings=case.get("n_assembly_rings"),
        n_core_rings=case.get("n_core_rings"),
    )
