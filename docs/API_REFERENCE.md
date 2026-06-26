# MOUSE API Reference

Function-level reference for the **MOUSE** (Microreactor Optimization Using
Simulation and Economics) routines, organized by subsystem. Each entry gives the
signature, what it reads/writes, the **underlying equation**, and **units**.

This document focuses on the routines that carry physics, geometry, engineering,
or economic equations. Pure OpenMC model-construction helpers and IO glue are
summarized rather than expanded.

---

## Conventions

**The `params` dict.** Most routines take a single mutable `params` dict and
*write their results back into it* (in place) rather than returning values. A
key written by one routine is an input to later ones. Entries below list the
keys a routine **Reads** and the keys it **Writes** (or the value it
**Returns**).

**Mass convention.** Densities are in **g/cm³** and geometry in **cm³**, so mass
in kilograms is almost always

```
mass [kg] = volume [cm^3] * density [g/cm^3] / 1000
```

**Unit cheat-sheet** (full glossary in Appendix A):

| Quantity | Unit |
|---|---|
| Length / radius / thickness / height | cm (engineering geometry); m inside heat-exchanger & pump correlations |
| Density | g/cm³ (neutronics/material tables); kg/m³ (thermo-hydraulic correlations) |
| Mass | kg |
| Thermal power | MWt; electric power MWe; mechanical/compressor power W (→ kW) |
| Temperature | K (absolute); °C where a correlation subtracts 273.15 |
| Enrichment, efficiency, capacity factor, rates, exponents | dimensionless fraction |
| Fuel lifetime, refueling, startup | days |
| Levelization period | years; **construction duration: months** |
| Burnup | MWd/kgHM |
| Cost | $ (escalated to `Escalation Year`); $/kW; LCOE $/MWh; LCOH $/MWth |

---

## Table of Contents

1. [Reactor Capability Registry](#1-reactor-capability-registry)
2. [Reactivity Control Elements](#2-reactivity-control-elements)
3. [Control Drums & Reflector/Moderator Masses](#3-control-drums--reflectormoderator-masses)
4. [Core Geometry Utilities](#4-core-geometry-utilities)
5. [Materials Database](#5-materials-database)
6. [Neutronics Post-Processing](#6-neutronics-post-processing)
7. [Vessel System](#7-vessel-system)
8. [Balance of Plant & Engineering Tools](#8-balance-of-plant--engineering-tools)
9. [Reactor Operation](#9-reactor-operation)
10. [Fuel Calculations](#10-fuel-calculations)
11. [Fuel-Lifetime Estimators](#11-fuel-lifetime-estimators)
12. [Parameter Builders](#12-parameter-builders)
13. [Estimate Service](#13-estimate-service)
14. [Cost Pipeline](#14-cost-pipeline)
15. [OpenMC Geometry Templates](#15-openmc-geometry-templates)
- [Appendix A: Units Glossary](#appendix-a-units-glossary)
- [Appendix B: Key Constants](#appendix-b-key-physical--calibration-constants)

---

## 1. Reactor Capability Registry

`core_design/reactor_registry.py` — describes each reactor concept by
independent *capability axes* instead of a bare name, so subsystems branch on
capabilities rather than on the reactor string.

### `class ReactorCapabilities` *(frozen dataclass)*
Capability descriptor for one concept. Fields (all dimensionless categorical):

| Field | Allowed values |
|---|---|
| `reactor_type` | `'LTMR'`, `'GCMR'`, `'HPMR'`, `'SRE'`, `'SFR'`, … |
| `spectrum` | `'thermal'`, `'fast'` |
| `moderator` | `None`, `'graphite'`, `'metal_hydride'` |
| `coolant_class` | `'gas'`, `'liquid_metal'`, `'heat_pipe'` |
| `control_element` | `'drum'`, `'rod'` |
| `vessel_architecture` | `'concentric_rvacs'`, `'pool'`, `'loop'` |
| `fuel_form` | `'triso'`, `'cylindrical_pin'`, `'metallic'` |
| `lattice` | `'hex'`, `'square'` |
| `pressurized` | `bool` (closed pressure vessel → full ellipsoidal heads) |
| `implemented` | `bool` (has a params builder in `reactor_config`) |

Predicates: `.is_fast`, `.is_moderated`, `.uses_drums`, `.uses_rods`.
Construction validates every axis against its allowed set.

### `get_capabilities(reactor_type) -> ReactorCapabilities`
Look up a concept; raises `KeyError` (listing known types) if absent.

### `is_registered(reactor_type) -> bool` / `implemented_reactor_types() -> tuple`
Membership test / the subset that `build_params` can actually build.

---

## 2. Reactivity Control Elements

`core_design/control_elements.py` — strategy interface that decouples
reactivity-control *sizing* from the drum-specific geometry, so rod-controlled
concepts plug in. Selected from the registry's `control_element` axis.

All implementations publish a common vocabulary:
`Control Element Type` (`'drum'`/`'rod'`), `Control Element Count`,
`Control Element Mass` [kg].

### `resolve_control_element(params) -> ControlElement`
Pick the strategy (explicit `params['Control Element Type']` → registry →
default `drum`), size it, and return it.

### `class DrumElement` → `.compute(params)`
Delegates to `drums.calculate_drums_volumes_and_masses` (§3) unchanged, then
mirrors `Control Drums Mass` onto the generic keys.

### `class RodElement` → `.compute(params)`
Control rods (absorber slug + clad tube + drive mechanism). **Reads:**
`Number of Control Rods` N (–), `Control Rod Radius` r_out [cm],
`Control Rod Clad Thickness` t_clad [cm], `Control Rod Absorber`/`Control Rod Clad`
(material names), `Active Height` [cm], optional `Control Rod Absorber Length`
L_abs [cm] (default = Active Height), `Control Rod Length` L_rod [cm] (default =
L_abs + 2·`Axial Reflector Thickness`), `Control Rod Drive Mass` m_drive [kg/rod].

**Equation** (r_abs = r_out − t_clad; ρ in g/cm³):
```
absorber_volume = pi * r_abs^2     * L_abs * N           [cm^3]
clad_volume     = pi * (r_out^2 - r_abs^2) * L_rod * N   [cm^3]
Control Rod Absorber Mass = absorber_volume * rho_abs  / 1000   [kg]
Control Rod Clad Mass     = clad_volume     * rho_clad / 1000   [kg]
Control Rod Drive Mass Total = N * m_drive                      [kg]
Control Element Mass = absorber + clad + drive_total            [kg]
All Rods Volume = pi * r_out^2 * L_rod * N   [cm^3]
All Rods Area   = pi * r_out^2 * N           [cm^2]
```
Validates N>0, 0 < t_clad < r_out, L_rod ≥ L_abs.

---

## 3. Control Drums & Reflector/Moderator Masses

`core_design/drums.py` — drum geometry/mass and the per-reactor reflector and
moderator mass models. Drum placement differs by reactor; the radius is
auto-maximized so drums just fit.

### `calculate_drums_volumes_and_masses(params)`
Resolve drum radius, then compute drum volumes and masses; for GCMR/HPMR also
auto-derive the dependent reflector/core geometry.

**Per-drum geometry** (r = `Drum Radius` [cm], t = `Drum Absorber Thickness`
[cm], h = `Drum Height` [cm]):
```
V_drum = pi * r^2 * h                                   [cm^3]
# absorber wedge (full annulus form):
V_absorber = pi * (r^2 - (r - t)^2) * h / 3             [cm^3]
# partial-coating form (if 'coating_angle' [deg] given):
V_absorber = [pi*r^2 - (pi/180)*angle*(r - t)^2] * h / 3
V_reflector = V_drum - V_absorber
```
**Totals & mass** (N = `Drum Count`; ρ from materials DB, g/cm³):
```
All Drums Volume = V_drum * N                           [cm^3]
All Drums Area   = All Drums Volume / h                 [cm^2]
Control Drum Absorber Mass  = V_absorber  * N * rho_absorber  / 1000  [kg]
Control Drum Reflector Mass = V_reflector * N * rho_reflector / 1000  [kg]
Control Drums Mass = absorber_mass + reflector_mass     [kg]
```

**Max-radius rules** (drum tube radius = r·(1+1/90) ≈ r·91/90 for LTMR/HPMR; r·46/45 for GCMR):
- **LTMR** `_calculate_max_ltmr_drum_radius` — drums on 6 hex faces, N ∈ {6,12,…,36}; bisection for the largest non-overlapping radius. Face positions: `r_face = apothem + drum_tube_radius`, drums spread along each edge.
- **GCMR** `_calculate_max_gcmr_drum_radius` — one drum per outer-ring hex cell: `r_max = (FTF/2)·(45/46)` [cm].
- **HPMR** `_calculate_max_hpmr_drum_radius` — drums on a ring of radius `r0 = (N_rings−1)·FTF + FTF/2`; no-overlap on the chord gives `max_tube = r0·sin(π/n)/(1−sin(π/n))`, `r_max = max_tube·(90/91)` [cm].

### `hexagonal_area_from_ftf(ftf) -> cm^2`
Regular hexagon area from flat-to-flat distance: `A = (√3/2)·ftf²`.

### `calculate_reflector_mass_LTMR / _GCMR(params)`, `calculate_reflector_and_moderator_mass_HPMR(params)`
Radial reflector = (circular core area − lattice/assembly area − `All Drums Area`)
× height × ρ/1000 [kg]; axial reflector = 2 × π·R_core²·t_axial × ρ/1000 [kg].
HPMR additionally computes moderator = (big-hex area − fuel area − heat-pipe
area) × H × ρ/1000 [kg]. (GCMR moderator handled separately below.)

### `calculate_moderator_mass_GCMR(params)`
Per-assembly moderator = hex area − fuel-compact area − coolant-channel area −
booster area, summed over all assemblies × H × ρ/1000 [kg]; booster mass summed
over annular booster shells. Fuel/coolant/booster counts follow the hex-lattice
combinatorics (`calculate_number_of_rings`).

### `calculate_moderator_mass(params) -> kg`
Simple pin moderator: `N_pins · π·r² · H · ρ / 1000`.

---

## 4. Core Geometry Utilities

`core_design/utils.py` — primitive geometry and hex-lattice combinatorics.
All radii/heights in **cm**, areas **cm²**, volumes **cm³**.

| Function | Equation |
|---|---|
| `circle_area(r)` | `π r²` [cm²] |
| `circle_perimeter(r)` | `2π r` [cm] |
| `cylinder_volume(r,h)` | `π r² h` [cm³] |
| `cylinder_radial_shell(r,h)` | `2π r h` [cm²] (lateral area) |
| `sphere_volume(r)` | `(4/3)π r³` [cm³] |
| `sphere_area(r)` | `4π r²` [cm²] |

### Hex lattice
- `calculate_hex_edge_length(params)` — `edge = pitch·(N_rings−1) + 0.6·pitch`,
  with `pitch = 2·FuelPinRadii[-1] + Pin Gap Distance` [cm].
- `calculate_hex_apothem(params)` — `apothem = sin(π/3)·edge = (√3/2)·edge` [cm].
- `calculate_core_radius_from_hex(params)` — `R_core = apothem + Radial Reflector Thickness` [cm].
- `calculate_number_of_rings(r)` — centered-hex count `= 3r² − 3r + 1` (r=1→1, 2→7, 3→19, 4→37). Equals GCMR `F_C`.
- `calculate_number_fuel_elements_hpmr(r)` — `N_total(r) − N_total(⌈r/2⌉)` (inner fuel, outer heat-pipe).
- `number_of_heatpipes_hmpr(params)` — writes per-assembly and total heat-pipe counts.

### Pin counting & particle counts
- `calculate_pins_in_assembly(params, pin_type)` — count of a pin label in the
  outer `N_rings` rows of `Pins Arrangement` (–).
- `calculate_total_number_of_TRISO_particles(params)` —
  `N_per_compact = ⌊PF · V_compact / V_particle⌋` with `V_compact = π·r_compact²·H`,
  `V_particle = (4/3)π·r_particle³`; total = per-compact × assemblies.

### Heat flux
- `calculate_heat_flux(params)` — pin surface flux:
  `q = P[MW] / (2π·r_pin·H·N_pins · 1e-4)` → **MW/m²** (1e-4 converts cm²→m²).
- `calculate_heat_flux_TRISO(params)` — kernel-surface flux:
  `q = P / (N·4π·r_kernel² · 1e-4)` → **MW/m²**.
- `monitor_heat_flux(params)` — warns if `Heat Flux` > `Heat Flux Criteria`.

### OpenMC orchestration (summary)
`run_openmc`, `run_depletion_analysis`, `openmc_depletion`,
`_run_isothermal_temperature_coefficients` drive 2-D OpenMC depletion and write
`keff 2D`, `keff 3D (2D corrected)`, leakage %, `Fuel Lifetime` [days],
`Mass U235/U238` [g], `Uranium Mass` [kg] (= (U235+U238)/1000). Reactivity
metrics: shutdown margin `SDM = (1 − k)/k · 1e5` [pcm]; isothermal temperature
coefficient `= Δ(1/k)/ΔT · 1e5` [pcm/K].

---

## 5. Materials Database

`core_design/openmc_materials_database.py`.

### `collect_materials_data(params) -> dict[str, openmc.Material]`
Builds every supported material (each in its own try/except, so missing params
skip a material rather than crash) and returns name → material, each carrying a
`.density` [g/cm³]. **Reads** `Enrichment`, `Common Temperature` [K], and
material-specific keys (`U_met_wo`, `H_Zr_ratio`, `UO2 atom fraction`, …).

Representative densities [g/cm³]: U_met 19.05, UO2 10.41, UC 13.0, UN 14.0,
UZr (U-10Zr) ~16.0, ZrH 5.6, YHx 4.28, Graphite 1.60, buffer_graphite 0.95,
PyC 1.9, SiC 3.18, B4C (nat/enr) 2.52, Be 1.84, BeO 3.01, SS304 7.93, NaK 0.85,
Helium 1.66e-4. Hydrides/graphite carry S(α,β) thermal-scattering tables.

---

## 6. Neutronics Post-Processing

### `core_design/peaking_factor.py :: compute_pin_peaking_factors(current_dir=".")`
From OpenMC depletion statepoints, tally pin kappa-fission power and compute the
power peaking factor per region and its envelope:
```
PF_i = P_i / mean(P)      [dimensionless]
Max_PF(step) = max_i PF_i
```
Returns `(summary_df[Step, Max_PF, Region], per_step_dict)`. Supports distribcell
(hex) and mesh (cylindrical) tallies.

### `core_design/correction_factor.py :: corrected_keff_2d(results, total_height, core_radius=None)`
Buckling-based axial/radial leakage correction of 2-D depletion k_eff, and fuel
cycle length. With 1-group collapsed cross-sections:
```
D   = 1 / (3 * Sigma_transport)            [cm]   (diffusion coefficient)
L^2 = D / Sigma_absorption                 [cm^2] (diffusion area)
B_axial^2  = (pi / (H_total + 2D))^2       [cm^-2]
B_radial^2 = (2.405 / (R_core + 2D))^2     [cm^-2]   (2.405 = first zero of J_0)
P_NL = 1 / (1 + L^2 * B^2)                 [dimensionless non-leakage prob.]
k_corrected = P_NL * k_2D
```
Leakage % = `(1 − P_NL)·100`. Fuel cycle length [days] = linear interpolation of
the corrected-k_eff curve to k=1. `total_height` = Active Height + 2·Axial
Reflector Thickness [cm].

---

## 7. Vessel System

`reactor_engineering_evaluation/vessel_stack.py` — generic ordered-layer vessel
model; `vessels_calcs.py` — `vessels_specs` entry point that selects an
architecture from the registry and writes the legacy mass/geometry keys.

### `class VesselLayer`
A cylindrical shell with an ellipsoidal bottom head: `name, radius [cm],
thickness [cm], bottom_depth [cm], material, head_factor, role`. Head factor:
`HEAD_FULL=1.0` (closed pressure vessel), `HEAD_HALF=0.5` (open-top pool),
`HEAD_NONE=0.0`.

### `compute_vessel_stack(layers, cylindrical_height, density_fn=materials_densities) -> VesselStackResult`
Integrate every layer with the **same** head+shell formula and sum (ρ in g/cm³):
```
head_vol = ellipsoid_shell(r, r, bottom_depth) * head_factor * thickness   [cm^3]
wall_vol = (circle_area(r + thickness) - circle_area(r)) * cyl_height       [cm^3]
layer_mass = (head_vol + wall_vol) * density(material) / 1000               [kg]
total_mass = sum(layer_mass)
outer_radius = last.radius + last.thickness
total_height = last.bottom_depth + cyl_height
```

### Architecture builders → `build_vessel_layers(params, architecture, pressurized)`
- `concentric_rvacs_layers` — 4 shells `inner → guard → cooling → intake`
  (RVACS air cooling); structural heads full when `pressurized` else half;
  reproduces the historical LTMR/GCMR/HPMR vessel masses exactly.
- `pool_layers` — `inner + guard` (single sodium pool + guard vessel), half heads.
- `loop_layers` — `inner` only (loop-type primary vessel).

### `vessels_specs(params)`
Compute the shared cylinder height, select architecture from
`get_capabilities(reactor type)`, build & integrate the stack, and write
`Vessel Mass`, `Guard/Cooling/Intake Vessel Mass`, `Total Vessels Mass`,
`Vessels Total Radius/Height`, `Vessel Height`, `Vessel Stack` (breakdown), and —
for a pressurized concentric stack — `RPV Outer Radius/Height`. All masses [kg],
lengths [cm].
```
Vessel Height = Active Height + 2*Axial Reflector Thickness
              + Vessel Lower Plenum Height + Vessel Upper Plenum Height
              + Vessel Upper Gas Gap                                  [cm]
```

### `ellipsoid_shell(a, b, c)` *(in tools.py; used as head surface area)*
Thomsen approximation of an ellipsoid surface area [cm²] (a,b,c in cm):
```
S = 4*pi * ( ((a*b)^1.6 + (a*c)^1.6 + (b*c)^1.6) / 3 )^(1/1.6)
```

---

## 8. Balance of Plant & Engineering Tools

`reactor_engineering_evaluation/tools.py` and `BOP.py`.

### `materials_densities(material) -> g/cm^3`
Structural-material lookup: stainless_steel/SS316 8.0, SS304 7.93,
low_alloy_steel/SA508 7.85, B4C 2.52, WEP 1.1.

### `material_specific_heat(material) -> J/(kg·K)`
Helium 5193, NaK 982, **Na/sodium 1270**.

### `cylinder_annulus_mass(r_out, r_in, h, material) -> kg`
`mass = 3.14·(r_out² − r_in²)·h · ρ/1000` (cm, g/cm³ → kg).

### `calculate_shielding_masses(params)`
In-vessel shield = annulus mass over `Vessel Height`; out-of-vessel shield outer
radius = `Vessels Total Radius` + thickness, mass × `Effective Density Factor`.
All [kg].

### `mass_flow_rate(params)`
Primary coolant mass flow from an energy balance:
```
m_dot = 1e6 * P[MWt] / (deltaT * cp)        [kg/s]
deltaT = Primary Loop Outlet T - Inlet T    [K]
cp = material_specific_heat(Coolant)        [J/(kg*K)]   (HPMR uses Secondary Coolant)
```
Writes `Coolant Mass Flow Rate` and `Primary Loop Mass Flow Rate` [kg/s];
divides by the per-loop load fraction when multiple loops are present.

### `compressor_power(params)` *(gas cycle)*
`P = ΔP[Pa]·ṁ[kg/s] / (η · ρ_He)` [W], with ρ_He = 3.3297 kg/m³ → `Primary Loop Compressor Power` [W].

### `compressor_wheel_diameter(params)` *(gas cycle)*
`D = (3.6/1.054) / (ΔP/ρ_He)^0.25 · √(V̇)` [m], `V̇ = ṁ/ρ_He` [m³/s].

### `GCMR_integrated_heat_transfer_vessel(params)`
PCHE volume (`Primary HX Mass / (ρ_HX·1e3) / 0.4`) + compressor cube → shell
volume × ρ → `Integrated Heat Transfer Vessel Mass` [kg] (0 if thickness 0).

### `BOP.calculate_heat_exchanger_mass(params) -> kg` *(printed-circuit HX / IHX)*
```
LMTD = (dT1 - dT2) / ln(dT1/dT2),  dT1=|Th_in - Tc_out|, dT2=|Th_out - Tc_in|   [K]
ht_area = P[MWt]*1e6 / (U * LMTD)            [m^2],  U = 500 W/m^2/K
n_channels = ht_area / channel_ht_area
alloy_volume = n_channels*pitch*thick - n_channels*(pi/8)*d_channel^2   [m^3]
hx_mass = alloy_volume * rho_ss              [kg],  rho_ss = 7850 kg/m^3
```
Reads `Primary/Secondary Loop Inlet/Outlet Temperature` [K].

### `BOP.calculate_primary_pump_mechanical_power(params)` *(liquid-metal loop)*
```
P_mech = m_dot * dP / (rho * eta)            [W]  -> /1000 -> kW
```
Defaults: `Primary Loop Pressure Drop` dP = 250 kPa, `Coolant Density` ρ = 750
kg/m³, `Pump Isentropic Efficiency` η = 0.75. Writes `Primary Pump Mechanical
Power` [kW].

### `BOP.calculate_secondary_pump_mechanical_power(m_dot) -> kW`
Static-head model: `P = ṁ·g·h/1000`, g = 9.81 m/s², h = 58.56·0.3048 m ≈ 17.85 m.

---

## 9. Reactor Operation

`reactor_engineering_evaluation/operation.py :: reactor_operation(params)`.
Converts lifetime/refueling/shutdown schedule into a capacity factor, annual
generation, and operator FTEs.
```
add_fuel_num = floor( 365 * Levelization Period / (Refueling Period + Fuel Lifetime) )
refuel_days/yr  = Refueling Period * add_fuel_num / Levelization Period
startup_days/yr = (after-refuel + after-emergency-shutdown contributions)
Capacity Factor = 1 - (refuel_days/yr + startup_after_refuel/yr + startup_after_shutdown/yr)/365
Annual Electricity Production = Capacity Factor * Power MWe * 365 * 24   [MWh]
FTE = days/yr * Work Hours Per Shift / Hours Per FTE
```
Times in **days**, period in **years**, generation in **MWh**.

---

## 10. Fuel Calculations

`reactor_engineering_evaluation/fuel_calcs.py`.

### `burnup_limited_fuel_lifetime_days(power_mwt, heavy_metal_mass_kg, discharge_burnup_mwd_per_kghm) -> days`
Generic burnup-limited lifetime (used by metallic-fuel concepts, e.g. SRE):
```
lifetime [days] = burnup [MWd/kgHM] * HM_mass [kg] / power [MWt]
```
Validates all inputs > 0.

### `fuel_calculations(params)`
Front-end enrichment material balance (natural-U feed, tails, SWU). **Reads**
`Mass U235`, `Mass U238` [g], `Enrichment` (–). With `U_mass = (U235+U238)/1000` [kg]:
```
nat_u_consum = U_mass * (E - 0.0025) / (0.0071 - 0.0025)   [kg]   (0.71% feed, 0.25% tails)
tail_waste   = nat_u_consum - U_mass                       [kg]
f(x) = (1 - 2x) * ln((1-x)/x)                               (value function)
SWU = U_mass*f(E) + tail_waste*f_tails - nat_u_consum*f_feed   [kg-SWU]
```
Writes `Natural Uranium Mass`, `Fuel Tail Waste Mass` [kg], `SWU`.

---

## 11. Fuel-Lifetime Estimators

`webapp/fuel_lifetime_estimator.py` (LTMR), `gcmr_fuel_lifetime_estimator.py`,
`hpmr_fuel_lifetime_estimator.py`. Each is a **K=4 distance-weighted KNN** over a
pre-computed OpenMC parametric study (XLSX), with features normalized to [0,1]
by training min/max and a physics-constrained normalized lifetime so off-grid
queries extrapolate sensibly. Lifetime in **days**; returns 0 for subcritical.

**LTMR** `estimate_ltmr_fuel_lifetime(n_rings_per_assembly, active_height, enrichment, power_mwt)`
```
N_pins = 3N^2 - 3N + 1
Lifetime = N_pins * H * A1 * (E - A2) / Power        [days]
  L* = Lifetime * Power / (N_pins * H)   (normalized; A1,A2 from OLS fit on K neighbors)
```
Features (E, N, H, P); training E∈[0.05,0.1975], N∈[6,24], H∈[50,180] cm, P∈[1,60] MWt.

**GCMR** `estimate_gcmr_fuel_lifetime(assembly_rings, core_rings, active_height, enrichment, power_mwt)`
```
F_A = 3(N_A-1)^2 - 3(N_A-1) + 1     (compacts per assembly)
F_C = 3 N_C^2 - 3 N_C + 1           (assemblies in core)
Lifetime = F_A * F_C * H * A1 * (E - A2) / Power     [days]
Uranium mass = 0.5776 * F_A * F_C * H   [g]   (0.5776 g per F_A*F_C*cm)
```
Features (E, N_A, N_C, H, P); training N_A∈[4,7], N_C∈[3,5], H∈[40,385] cm.

**HPMR** `estimate_hpmr_fuel_lifetime(n_rings_per_assembly, n_rings_per_core, active_height, enrichment, power_mwt)`
```
N_pins(N_A,N_C) = (3 N_C^2 - 3 N_C) * f_HPMR(N_A),  f_HPMR(r) = N_total(r) - N_total(ceil(r/2))
L* = LT * P / (N_pins * H * E)                 (normalized; ~constant)
hpmr_total_uranium_mass_g = 1.6116 * N_pins * H   [g]   (1.6116 g per pin per cm)
```
Features (E, N_C, H, P), N_A fixed at 6; training H∈[136,1056] cm, P∈{1,5,20,60} MWt.

**Shared helpers** (all three): `_..._knn_scalar(column, …)` distance-weighted
interpolation of any parametric column; `get_..._peaking_factor` (–),
`get_..._axial/total_leakage_pct` (%), `get_..._keff_curve(…, anchor_lifetime_days)`
(times [days], k_eff [–]). **Physics leakage fallback** (out-of-range geometry):
```
delta = savings * reflector_thickness          (savings: LTMR 0.55, GCMR/HPMR 0.65)
B_radial^2 = (2.405/(R+delta_r))^2,  B_axial^2 = (pi/(H+2*delta_z))^2
P_NL = 1 / (1 + M^2 * B_total^2)               (M^2: LTMR 60, GCMR/HPMR 220 cm^2)
total_leak% = (1-P_NL)*100,  axial_leak% = total% * B_axial^2/B_total^2
```

---

## 12. Parameter Builders

`webapp/reactor_config.py` — assembles a fully-populated `params` dict for a
reactor (no OpenMC), ready for the cost engine.

### `build_params(reactor_type, power_mwt, enrichment, user_overrides, n_rings_per_assembly=None, active_height=None, n_assembly_rings=None, n_core_rings=None) -> dict`
Validate the per-reactor required inputs, dispatch to `_build_<type>`, then apply
`user_overrides` (which win). Raises `NotImplementedError` for a registered but
unbuilt concept (e.g. SFR) and `ValueError` for an unknown one. **Units:**
`power_mwt` MWt, `enrichment` fraction, `active_height` cm, ring counts –.

### `_build_ltmr / _build_gcmr / _build_hpmr / _build_sre (params)`
Populate ~10 sections (materials, geometry, control, overall system, fuel
lifetime, primary loop/BoP, shielding, vessels, operation, buildings/economics).
LTMR/GCMR/HPMR use control **drums** + concentric RVACS vessels + KNN lifetime;
**`_build_sre`** uses control **rods** (via `resolve_control_element`) + **loop**
vessel + sodium + metallic fuel + burnup-limited lifetime:
```
Uranium Mass = N_rods * pi * r_meat^2 * H * Fuel Density / 1000   [kg]
Fuel Lifetime = burnup_limited_fuel_lifetime_days(P, Uranium Mass, Discharge Burnup)
Moderator/Reflector Mass = graphite_volume * rho_graphite / 1000  [kg]
```

### `interpolate_openmc_results(reactor_type, power_mwt, enrichment) -> dict`
Bilinear interpolation of a parametric-study CSV: `Mass U235/U238` [g]
(1-D in enrichment), `Fuel Lifetime` [days] (enrichment × power); raises
`SubcriticalError` if 0.

---

## 13. Estimate Service

`webapp/estimate_service.py` — app-facing orchestration (no Streamlit).

### `EstimateInputs` *(frozen dataclass)*
One design point. Key fields/units: `power_mwt` MWt, `enrichment` fraction,
`interest_rate`/`discount_rate` fraction/yr, **`construction_duration` months**,
`debt_to_equity` ratio, `emergency_shutdowns` /yr, `startup_duration(_refueling)`
days, `tax_credit_type` ITC/PTC, `tax_credit_value` fraction (ITC) or $/MWh (PTC),
`plant_lifetime` years, `active_height` cm, ring counts –.

### `run_estimate(inputs) -> EstimateResult`
`build_params` → `bottom_up_cost_estimate('cost/Cost_Database.xlsx', params)` →
`cost_drivers_estimate` → `transform_dataframe`. Returns `display_df`,
`enriched_df` (with `FOAK/NOAK LCOE` [$/MWh]), `detailed_sorted_df`, `params`.

### `run_lcoe_at_noak_unit(inputs) -> LcoeAtNoakResult`
LCOE anchor at a given `noak_unit_number`; returns `mean`, `std` [$/MWh] and
optional diagnostics. `_get_mean_std(df, account, which)` pulls the mean/std cost
[$] for an account from the engine output.

---

## 14. Cost Pipeline

`cost/` — bottom-up code-of-account model. Costs are escalated to
`Escalation Year` and reported in **$**, `$/kW`, or **$/MWh**.

### Escalation — `cost/cost_escalation.py`
- `calculate_inflation_multiplier(file, base_year, cost_type, escalation_year)` —
  `multiplier = index(escalation_year)/index(base_year)` (1 if type `'NA'`).
- `escalate_cost_database(file, escalation_year, params, sheet_name)` —
  resolve param-referenced costs, multiply by the inflation factor
  (`Adjusted = Original × multiplier`), and merge the *Economics Parameters*
  sheet into `params`.

### Account scaling — `cost/cost_scaling.py`
- `scale_cost(db, params)` — per row, with X = `params[Scaling Variable]`,
  X_ref = `Scaling Variable Ref Value`, exponent n:
  ```
  standard (X_ref > 0):  cost = fixed + unit * X^n / X_ref^(n-1)
  standard (X_ref <= 0): cost = fixed + unit * X            (linear)
  X = 0 (and X is a key):  cost = 0
  nonstandard:           cost = non_standard_cost_scale(...)
  ```
  With `Number of Samples > 1`, `fixed`/`unit` are sampled (Lognormal/Uniform)
  and `n` from a Truncated Normal.
- `non_standard_cost_scale(account, unit_cost, X, n, params)` — per-account
  engineering correlations, e.g.:
  ```
  pumps 222.11/12:  cost = [0.2/(1 - Pump Isentropic Efficiency) + 1] * unit * X^n
  compressor 222.13: cost = [((T_out[K]-273.15)/650)^1.29 * (P_comp[W]/1e6/2.6)^0.74] * unit
  enrichment 253:   cost = premium * unit * X^n   (premium 1.0 if E<0.10, 1.15 if 0.10<=E<0.20)
  staffing 711/712/713/81: cost = FTE_multiplier * unit * X^(±n)
  ```
- `scale_redundant_BOP_and_primary_loop(df, params)` — multiply 222.x by
  `Primary Loop Count`, 232.x/213.1 by `BoP Count`, 226 by `Primary Loop Purification`.

### FOAK→NOAK learning — `cost/cost_estimation.py`
- `learning_rate_multiplier(learning_rate, n_units)` — **Wright's law**:
  `multiplier = (1 − learning_rate)^log2(min(n_units, 100))` (capped at the 100th unit).
- `FOAK_to_NOAK(df, params)` — `NOAK = FOAK × multiplier` per `FOAK to NOAK
  Multiplier Type` (No Learning, Licensing, Factory Primary Structure, Factory
  Drums, …; onsite-learning units = 2 × `NOAK Unit Number`).
- `bottom_up_cost_estimate(file, params)` — full pipeline: validate tax credits →
  escalate → `remove_irrelevant_account` → `reactor_operation` → per sample
  {scale → redundant-BoP → FOAK→NOAK → roll up base (acct 1–2) → indirect
  31/32/75/82 → decommissioning 78 → roll up 3–5 → interest 62 → roll up 6 → TCI
  → roll up 7–8 → LCOE} → mean & std across samples. Returns the account table
  ($, with std columns).
- Roll-up (`update_high_level_costs` / `calculate_high_level_accounts_cost`):
  parent cost = Σ children, leaf-to-root over levels 4→0.

### Indirect, financing, decommissioning, LCOE — `cost/non_direct_cost.py`
- `_crf(rate, period)` — **capital recovery factor**:
  `CRF = rate·(1+rate)^period / ((1+rate)^period − 1)` (= 1/period if rate=0). [1/yr]
- `calculate_accounts_31_32_75_82_cost` — indirect field/owner costs (ratios of
  direct cost); annualized component **replacement** = capital × `CRF(Discount
  Rate, replacement_period_yr)`; annual fuel (acct 82) = fuel × `CRF(Discount
  Rate, refueling_period_yr)`. Replacement/refueling periods built from
  `Fuel Lifetime` + `Refueling Period` + startup [days] ÷ 365.
- `calculate_high_level_capital_costs` — `OCC = Σ(accts 10,20,30,40,50)` [$];
  `OCC per kW = OCC/(1000·MWe)`; `OCC excl. fuel = OCC − acct 25`. **Interest
  during construction** (acct 62), with debt fraction `d = D:E/(1+D:E)` and
  construction duration in **months**:
  ```
  B = 1 + exp( ln(1+rate) * months/12 )
  C = ( ln(1+rate) * (months/12) / pi )^2 + 1
  Interest = d * OCC * (0.5 * B / C - 1)      [$]
  ```
- `calculate_decommissioning_cost` (acct 78) — future value =
  `(acct10+acct20) × CAPEX-to-Decommissioning Ratio` (default 0.15), annualized
  over `Levelization Period` at `Annual Return`.
- `calculate_TCI` — `TCI = OCC + acct 60` [$]; `TCI per kW = TCI/(1000·MWe)`;
  optional ITC-adjusted OCC/TCI when `ITC credit level` is set (and the unit is
  within the claiming-units cutoff).
- `energy_cost_levelized(params, df)` — **LCOE** by discounted cash flow over
  `Levelization Period` years at `Discount Rate` d:
  ```
  LCOE = [ TCI(yr0) + sum_{i=1..LP} (O&M_i + fuel_i)/(1+d)^i ]
         / [ sum_{i=1..LP} E_i/(1+d)^i ]            [$/MWh]
  E_i = Power MWe * Capacity Factor * 365 * 24      [MWh]
  ```
  PTC: subtract `PV(PTC credit)/PV(energy)` (grossed up by `1/(1−Tax Rate)`).
  LCOH: scaled capital/O&M × `Thermal Efficiency` → `$/MWth`.

### Account utilities — `cost/code_of_account_processing.py`
- `remove_irrelevant_account(df, params)` — drop a row unless its
  `Optional Variable` is present in `params` and equals (or is contained in)
  its `Optional Value` (same for `Sec Optional Variable`). This is how
  reactor-type / material / coolant gating selects applicable accounts.
- `find_children_accounts(df)` — attach each parent's child row indices.
- `get_estimated_cost_column(df, 'F'|'N'|'F std'|'N std')` — locate a cost column.
- `create_cost_dictionary(df, params, tracked)` — flatten tracked params + OCC/
  TCI/LCOE (+ ITC/PTC variants) for CSV export.

### Sampling — `cost/sampling.py`
`sampler(distribution, **kwargs)` dispatches:
- Lognormal — μ, σ from ln(low, class3, high); returns one draw.
- Truncated Normal — rejection-sample N(mean, std) into [lower, upper].
- Uniform — U(low, high).

### Cost drivers / standalone LCOE — `cost/cost_drivers.py`, `cost_drivers_and_markets/lcoe.py`
`energy_cost_levelized_per_acct(params, capital_cost, ann_cost)` applies the same
DCF to one account → its LCOE contribution [$/MWh]; `cost_drivers_estimate`
enriches the table with `FOAK/NOAK LCOE` columns and optional bar charts.
`cost_drivers_and_markets/lcoe.py :: energy_cost_levelized(...)` is the standalone
(no tax-credit) DCF LCOE used by legacy analyses.

---

## 15. OpenMC Geometry Templates

`core_design/openmc_template_{LTMR,GCMR,HPMR}.py`,
`openmc_template__HPMR_vtb.py`, `pins_arrangement.py` — build the OpenMC
neutronics model (geometry, materials, tallies, settings) for each reactor.
**These require OpenMC installed** (unlike the params/cost path). Each
`build_openmc_model_*(params)` writes `materials.xml`, `geometry.xml`,
`tallies.xml`, `settings.xml`. Lengths in **cm**, angles in **degrees**.

**Geometry hierarchy** (all hex-based): pin cells → hex lattice (rings per
assembly) → assemblies → hex lattice (core rings) → radial reflector + control
drums → outer boundary (vacuum); axially: bottom reflector | active height | top
reflector.

### `build_openmc_model_LTMR(params)` *(openmc_template_LTMR.py:433)*
Cylindrical fuel/moderator pins in hex assemblies; drums on the 6 hex faces.
- **Reads:** `Fuel Pin Radii` [cm] (cumulative: insert, gap, fuel meat, gap,
  clad), `Moderator Pin Radii` [cm], `Pin Gap Distance` [cm], `Pins Arrangement`,
  `Number of Rings per Assembly`, `Number of Drums` (6/12/…/36), `Drum Radius`
  [cm] (auto-max if absent), `Drum Absorber Thickness` [cm], `Drum Absorber Arc
  Degrees` [deg], `Active Height` [cm], material names, `Common Temperature` [K],
  `cross_sections_xml_location`; flags `Shutdown Margin Calc`, `Isothermal
  Temperature Coefficients`, `plotting`.
- **Key formulas:** pin pitch `= 2·FuelPinRadii[-1] + Pin Gap Distance`; apothem
  `= Assembly FTF/√3`; drum tube radius `= Drum Radius·(1 + 1/90)`; drum centers
  at radius `apothem + drum_tube_radius` on each of the 6 face normals.

### `build_openmc_model_GCMR(params)` *(openmc_template_GCMR.py:13)*
TRISO particles stochastically packed (`openmc.model.pack_spheres`) in compacts;
6 drums in the outer-ring hex cells.
- **Reads:** `Fuel Pin Radii` [cm] (TRISO shells: kernel, buffer, PyC, SiC, PyC),
  `Compact Fuel Radius` [cm], `Packing Fraction` (–), `Matrix Material`,
  `Lattice Pitch` [cm], `Assembly Rings`, `Core Rings`, `Coolant Channel Radius`
  [cm], `Moderator Booster Radii` [cm], `Drum Radius`/`Drum Absorber Thickness`
  [cm], `Active Height` [cm], …
- **Key formulas:** hex lattice radius `= Lattice Pitch/√3`; absorber arc 120°
  (plane coefficient `1/√3`); fuel-kernel volume `(4/3)·π·FuelPinRadii[0]³`.

### `build_openmc_model_HPMR(params)` *(openmc_template_HPMR.py:400)*
Fuel pins + heat pipes alternating in hex assemblies; drums on a placement ring.
- **Reads:** `Fuel Pin Radii` / `Heat Pipe Radii` [cm], `Lattice Pitch` [cm],
  `Number of Rings per Assembly` / `per Core`, `Assembly FTF` [cm], `hexagonal
  Core Edge Length` [cm], `Drum Count` (6/12/18/24), `Drum Radius`/`Drum Absorber
  Thickness` [cm], `Active Height` [cm], …
- **Key formulas:** drum ring radius `= (N_core_rings−1)·FTF + FTF/2 +
  drum_tube_radius`; angular spacing `360°/Drum Count`; absorber snapped to the
  nearest hex-face normal (nearest 60°), +180° in shutdown-margin mode.

### `update_ltmr_reflector_geometry_from_drums(params)` *(openmc_template_LTMR.py)*
Derive `Core Radius`, `Radial Reflector Thickness`, `Axial Reflector Thickness`
[cm] from the drum layout (drums drive reflector size). Called by the LTMR param
builder before the drum-mass calculation.

### Energy groups (all templates)
Shared **11-group** MGXS structure [eV] (LTMR:626, GCMR:440, HPMR:511):
```
[1e-5, 6.7e-2, 3.2e-1, 1, 4, 9.88, 4.81e1, 4.54e2, 4.9e4, 1.83e5, 8.21e5, 4e7]
```
Thermal-optimized (HPMR report, table 5); a fast-spectrum concept would need a
different group structure.

### Control drum construction & placement
| Template | Count | Placement | Absorber arc | Clearance gap |
|---|---|---|---|---|
| LTMR | 6–36 | 6 hex faces, N/6 per face | parametric [deg] (two cut planes) | r/90 |
| GCMR | 6 | outer-ring hex cells (60° symmetry) | 120° (`1/√3`) | r/45 |
| HPMR | 6–24 | placement ring, snapped to 60° | per-drum plane | r/90 |
| HPMR-vtb | 12 | 6 hex vertices | parametric `coating_angle` | inner/outer shells |

### Tallies
11-group MGXS (absorption, diffusion coefficient, transport, scatter matrix,
total) + pin power as **`kappa-fission`**: LTMR/HPMR via `DistribcellFilter` on
the fuel-meat cell; GCMR via a 20×20 `MeshFilter` + `MaterialFilter` on the
kernel.

### `pins_arrangement.py :: LTMR_pins_arrangement`
A nested list — one inner list per hex ring, each `[...] * 6` for the six
sectors — of `'FUEL'`/`'MODERATOR'` placeholders defining the LTMR pin layout.
The builder takes the outer `Number of Rings per Assembly` rings and maps
`'FUEL'`→fuel-pin universe, `'MODERATOR'`→moderator-pin universe for
`openmc.HexLattice`. (No `'CONTROL_ROD'` placeholder — rod-controlled concepts
would add one.)

### `openmc_template__HPMR_vtb.py :: class OpenMC_HPMR`
Detailed HPMR test-bed generator. `gen_input(parms, Tdict, rundir, …)` builds a
parametric model with axially-divided depletable fuel and 12 staggered drums;
`run_nominal(…)` runs drum orientations 0°/90°/180° (the last with depletion);
`postproc()` extracts shutdown margin, lifetime (burnup to k=1), peaking factors
(Fq, Fdh), and linear heat rate to CSV.

---

## Appendix A: Units Glossary

| Symbol / key | Meaning | Unit |
|---|---|---|
| `Power MWt` / `Power MWe` | thermal / electric power | MWt / MWe |
| `Enrichment` | U-235 atom fraction | – (0–1) |
| radii, thicknesses, heights, `*FTF`, apothem | geometry | cm |
| areas / volumes | geometry | cm² / cm³ |
| density (`.density`, `materials_densities`) | material | g/cm³ |
| density (`Coolant Density`, ρ_ss) | thermo-hydraulic | kg/m³ |
| masses (`* Mass`) | component mass | kg |
| `Common Temperature`, loop temperatures | absolute temperature | K |
| `cp` (specific heat) | thermal | J/(kg·K) |
| `* Mass Flow Rate` | coolant flow | kg/s |
| `* Mechanical Power`, `Compressor Power` | pump/compressor | W or kW |
| `Heat Flux` | surface flux | MW/m² |
| `Fuel Lifetime`, `Refueling Period`, startup | time | days |
| `Levelization Period` | plant life | years |
| `Construction Duration` | build time | **months** |
| `Discharge Burnup` | fuel burnup | MWd/kgHM |
| `Capacity Factor`, efficiencies, rates, exponents | – | dimensionless |
| `Annual Electricity Production` | generation | MWh |
| costs (`* Cost`, OCC, TCI) | money | $ (at `Escalation Year`) |
| `* per kW` | specific capital | $/kW |
| LCOE / LCOH | levelized cost | $/MWh / $/MWth |
| reactivity (SDM, temp coeff) | – | pcm, pcm/K |

## Appendix B: Key Physical & Calibration Constants

| Constant | Value | Where |
|---|---|---|
| KNN neighbors `K` | 4 | all fuel-lifetime estimators |
| Migration area M² (LTMR / GCMR / HPMR) | 60 / 220 / 220 cm² | physics leakage fallback |
| Reflector savings (LTMR / GCMR-HPMR) | 0.55 / 0.65 | physics leakage fallback |
| Uranium per pin (HPMR) | 1.6116 g/(pin·cm) | `hpmr_total_uranium_mass_g` |
| Uranium per F_A·F_C·cm (GCMR) | 0.5776 g | GCMR uranium mass |
| Drum tube clearance (LTMR/HPMR / GCMR) | ×91/90 / ×46/45 | drum sizing |
| HX overall U | 500 W/m²/K | `calculate_heat_exchanger_mass` |
| ρ stainless (HX) | 7850 kg/m³ | `calculate_heat_exchanger_mass` |
| Pump defaults | ΔP 250 kPa, ρ 750 kg/m³, η 0.75 | `calculate_primary_pump_mechanical_power` |
| He density (gas cycle) | 3.3297 kg/m³ | `compressor_power` |
| Wright's-law cap | 100th unit | `learning_rate_multiplier` |
| Decommissioning ratio (default) | 0.15 | `calculate_decommissioning_cost` |
| Radial buckling root (J₀) | 2.405 | leakage / `corrected_keff_2d` |

---

*Generated as a developer reference for MOUSE. Equations are transcribed from the
source; see the cited modules for exact implementation and assumptions.*
