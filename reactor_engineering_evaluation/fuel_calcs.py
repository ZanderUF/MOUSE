# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED

import numpy as np


def burnup_limited_fuel_lifetime_days(power_mwt, heavy_metal_mass_kg,
                                      discharge_burnup_mwd_per_kghm):
    """Full-power days until the fuel reaches its discharge-burnup limit.

        lifetime [days] = burnup [MWd/kgHM] * HM mass [kg] / power [MWt]

    This is the natural fuel-lifetime model for metallic-fuel concepts (e.g. an
    SRE-like or sodium fast reactor) that have no surrogate burnup-vs-geometry
    training data: the achievable discharge burnup is supplied as a design input
    (for unalloyed/low-alloy metallic uranium it is typically swelling-limited).
    """
    if power_mwt <= 0:
        raise ValueError(f"power_mwt must be positive, got {power_mwt}.")
    if heavy_metal_mass_kg <= 0:
        raise ValueError(
            f"heavy_metal_mass_kg must be positive, got {heavy_metal_mass_kg}.")
    if discharge_burnup_mwd_per_kghm <= 0:
        raise ValueError(
            "discharge_burnup_mwd_per_kghm must be positive, got "
            f"{discharge_burnup_mwd_per_kghm}.")
    return discharge_burnup_mwd_per_kghm * heavy_metal_mass_kg / power_mwt


def fuel_calculations(params):
    
    """Front-end enrichment material balance (natural-U feed, tails, separative work).

    With U_mass = (Mass U235 + Mass U238)/1000 [kg] and enrichment E (–):
      nat_u_consum = U_mass*(E - 0.0025)/(0.0071 - 0.0025)   [kg]  (0.71% feed, 0.25% tails)
      tail_waste   = nat_u_consum - U_mass                    [kg]
      value fn f(x) = (1 - 2x)*ln((1-x)/x)
      SWU = U_mass*f(E) + tail_waste*5.96 - nat_u_consum*4.87   [kg-SWU]
    Mutates params: writes 'Natural Uranium Mass', 'Fuel Tail Waste Mass' [kg], 'SWU'.
    """
    U_mass = (params['Mass U235'] + params['Mass U238']) / 1000  # mass of uranium only (kg)
    nat_u_consum = U_mass*(params['Enrichment'] -0.0025)/(0.0071-0.0025)  # kg
    tail_waste = nat_u_consum - U_mass  # kg

    # Value functions
    f_val_fun = (1-2*params['Enrichment'])*np.log((1-params['Enrichment'])/params['Enrichment'])
    tail_waste_val_fun = 5.96
    nat_u_waste_val_fun = 4.87


    kg_SWU = (U_mass*f_val_fun+tail_waste*tail_waste_val_fun- nat_u_consum *nat_u_waste_val_fun)
    
    params['Natural Uranium Mass'] = nat_u_consum
    params['Fuel Tail Waste Mass'] = tail_waste  # kg
    params['SWU'] = kg_SWU
 