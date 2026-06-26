# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED
import numpy as np 

def ellipsoid_shell(a, b, c):
    """Thomsen approximation of an ellipsoid surface area (ellipsoidal vessel-head area).
    S = 4*pi*( ((a*b)^1.6 + (a*c)^1.6 + (b*c)^1.6)/3 )^(1/1.6)
    Inputs: a, b, c semi-axes in cm. Returns surface area in cm^2.
    """
    return 4*np.pi*np.power(((a*b)**1.6 + (a*c)**1.6 + (b*c)**1.6)/3, 1/1.6)

def circle_area(r):
    """Area of a circle.
    area = pi*r^2
    Input: r radius in cm. Returns area in cm^2.
    """
    return (np.pi) * r **2


def materials_densities(material):
    """Structural-material density lookup by name.
    Returns density in g/cm^3 (stainless_steel/SS316 8.0, SS304 7.93,
    low_alloy_steel/SA508 7.85, B4C 2.52, WEP 1.1).
    Input: material name (str). Returns density in g/cm^3.
    """
    material_densities = {
    "stainless_steel": 8.0,  # Approximate density of stainless steel
    "SS316": 8.0,            # Approximate density of SS316
    "SS304": 7.93,           # Approximate density of SS304
    "low_alloy_steel": 7.85, # Approximate density of SA508 Gr3 Cls 1
    "SA508": 7.85,           # Approximate density of SA508 Gr3 Cls 1
    "B4C_enriched": 2.52,    # Approximate density of boron carbide
    "B4C_natural": 2.52,     # Approximate density of boron carbide
    "WEP": 1.1,              # WEP density (water extended polymer)
    }
    return material_densities[material] # in gram/cm^3

def material_specific_heat(material):
    """Coolant specific-heat lookup by name.
    Returns specific heat in J/(kg*K) (Helium 5193, NaK 982, Na/sodium 1270).
    Input: material name (str). Returns specific heat in J/(kg*K).
    """
    material_cp = {
        "Helium": 5193,  # J/(kg·K)
        "NaK": 982.,     # J/(kg·K)
        "Na": 1270.,     # J/(kg·K) — liquid sodium at ~450 °C (SRE/SFR primary coolant)
        "sodium": 1270., # alias for Na
    }
    return material_cp[material]  # J/(kg·K)

def cylinder_annulus_mass(outer_radius , inner_radius,height, material ):
    """Mass of an annular (hollow) cylinder for a given material.
    volume = 3.14*(outer_radius^2 - inner_radius^2)*height [cm^3]; mass = volume*density/1000 [kg]
    Inputs: outer_radius, inner_radius, height in cm; material name (str). Returns mass in kg.
    """
    volume = 3.14* (outer_radius**2 - inner_radius**2) * height
    mass = volume* materials_densities(material)/1000  # kg
    return mass # in kg

def calculate_shielding_masses(params):
    """Compute in-vessel and out-of-vessel shield masses and write them into params.
    Uses cylinder_annulus_mass over the vessel height; the outer-shield mass is scaled
    by 'Out Of Vessel Shield Effective Density Factor'.
    Mutates params: writes 'In Vessel Shield Mass' and 'Out Of Vessel Shield Mass' [kg].
    """
    params['In Vessel Shield Mass'] = cylinder_annulus_mass(params['In Vessel Shield Outer Radius'],\
    params['In Vessel Shield Inner Radius'], params['Vessel Height'], params['In Vessel Shield Material'] )
    params['Outer Shield Outer Radius'] = params['Out Of Vessel Shield Thickness'] + params['Vessels Total Radius']
    params['Outer Shield Inner Radius'] = params['Outer Shield Outer Radius'] - params['Out Of Vessel Shield Thickness']

    outer_shield_mass = cylinder_annulus_mass(params['Outer Shield Outer Radius'], params['Outer Shield Inner Radius'],\
    params['Vessels Total Height'], params['Out Of Vessel Shield Material']) 
    params['Out Of Vessel Shield Mass'] = params['Out Of Vessel Shield Effective Density Factor'] * outer_shield_mass

def mass_flow_rate(params):
    """Primary-coolant mass flow rate from a steady-state energy balance.
    m_dot = 1e6*Power_MWt/(deltaT*cp) [kg/s], deltaT = Outlet - Inlet Temperature [K],
    cp = material_specific_heat(Coolant) [J/(kg*K)] (HPMR uses Secondary Coolant).
    Mutates params: writes 'Coolant Mass Flow Rate' and 'Primary Loop Mass Flow Rate' [kg/s].
    """
    loop_factor = 1
    thermal_power_MW = params['Power MWt']
    if 'Primary Loop per loop load fraction' in params.keys():
        loop_factor = params['Primary Loop per loop load fraction']
        thermal_power_MW = params['Power MWt'] * loop_factor
        
    deltaT =  params['Primary Loop Outlet Temperature'] - params['Primary Loop Inlet Temperature']
    if params['reactor type'] == "HPMR":
        coolant = params['Secondary Coolant']
    else:    
        coolant = params['Coolant']
    coolant_specific_heat = material_specific_heat(coolant)
    m_dot = 1e6 * thermal_power_MW/ (deltaT * coolant_specific_heat)
    params['Coolant Mass Flow Rate']  = m_dot / loop_factor # For Reactor Mass Flow Rate
    params['Primary Loop Mass Flow Rate'] = m_dot # For individual Primary Loop Mass Flow Rate
    
def compressor_power(params):
    """Gas-cycle compressor power from pressure drop and isentropic efficiency.
    power = Primary Loop Pressure Drop [Pa] * Primary Loop Mass Flow Rate [kg/s]
            / (Compressor Isentropic Efficiency * rho_he), rho_he = 3.3297 kg/m^3.
    Mutates params: writes 'Primary Loop Compressor Power' [W].
    """
    # Estimates the required compressor power based on a simplified
    # model using pressure drop and compressor isentropic efficiency

    rho_he = 3.3297  # kg/m3. TODO: Consider importing CoolProp to estimate density based on cold-leg temperature and pressure
    power = params['Primary Loop Pressure Drop']*params['Primary Loop Mass Flow Rate']/params['Compressor Isentropic Efficiency']/rho_he
    params['Primary Loop Compressor Power'] = power # W
    return

def compressor_wheel_diameter(params):
    """Approximate compressor wheel diameter from specific-diameter scaling.
    diameter = (3.6/1.054)/(dP/rho_He)^0.25 * sqrt(Vdot), Vdot = mdot/rho_He [m^3/s],
    rho_He = 3.330 kg/m^3, dP = Primary Loop Pressure Drop [Pa]. Returns diameter in m.
    """
    # Estimates the approximate compressor size based on its specific
    # diameter, matched to the MIGHTR horizontal HTGR design.
    # Ref for specific diameter:
    #  https://www.dropbox.com/scl/fi/fnqdg2hyi6y4ozu9p7nyu/final-report-str-mech-ARDP-redacted-V3.pdf?rlkey=h97dii28tvf0bxtffo8q62tn5&st=zsls1bs2&dl=0
    ref_specific_diameter = 3.6  # dimensionless
    rho_He = 3.330  # kg/m3 for He at 4 MPa, 300 °C. TODO: use a He density correlation or CoolProp to estimate density based on cold-leg temperature and pressure
    Vdot_gcmr = params['Primary Loop Mass Flow Rate'] / rho_He  # m3/s — volumetric flow rate
    dP = params['Primary Loop Pressure Drop']
    diameter = ref_specific_diameter/1.054 / (dP/rho_He)**0.25 * np.sqrt(Vdot_gcmr) # m
    return diameter

def GCMR_integrated_heat_transfer_vessel(params):
    """Size the integrated heat-transfer vessel housing the PCHE and circulator.
    PCHE volume = Primary HX Mass/(rho_HX*1e3)/0.4; compressor volume ~ wheel_diameter^3.
    Mutates params: writes 'Integrated Heat Transfer Vessel Mass' [kg] and Outer Volume [m^3]
    (both set to 0 when Integrated Heat Transfer Vessel Thickness is 0).
    """
    # Calculates the required parameters for the
    # GCMR Integrated Heat Transfer Vessel that houses:
    #   circulator, PCHE, piping, valves, insulation

    contingency = 0.3  # accounts for the volume/mass of valves, fittings, and connections
    PCHE_volume = (params['Primary HX Mass'] / (materials_densities(params['HX Material'])*1e3) / 0.4)  # accounts for assumed 60% coolant channel void fraction
    compressor_volume = (compressor_wheel_diameter(params))**3  # approximated as a cube with side = wheel diameter. TODO: improve compressor sizing

    vessel_inner_volume = (1+contingency)*(PCHE_volume + compressor_volume)  # assumes a cube-like structure
    vessel_outer_volume = (vessel_inner_volume**(1/3)+ 1e-2*params['Integrated Heat Transfer Vessel Thickness'])**3  # m3
    vessel_volume = vessel_outer_volume - vessel_inner_volume
    vessel_density = materials_densities(params['Integrated Heat Transfer Vessel Material'])*1e3
    params['Integrated Heat Transfer Vessel Outer Volume'] = vessel_outer_volume
    params['Integrated Heat Transfer Vessel Mass'] = vessel_volume * vessel_density
    
    if params['Integrated Heat Transfer Vessel Thickness'] == 0:
        params['Integrated Heat Transfer Vessel Outer Volume'] = 0
        params['Integrated Heat Transfer Vessel Mass'] = 0

    # Rough estimate of the mass supported by the support structure:
    # Primary HX + Integrated Heat Transfer Vessel + Compressor + Valves/Fittings/Bolts/etc.
    # params['Integrated Heat Transfer System Mass'] = params['Primary HX Mass'] + (vessel_volume * vessel_density) + compressor_volume*8000