# -*- coding: utf-8 -*-
"""
Classes to contain RocketMotor objects

Created on Sat Feb 15 13:30:55 2025

@author: hugo
"""

from numpy import where

#%% Constants ----------------------------------------------------------------
g = 9.81        # standard gravity acceleration
R_prime = 8.3144598 # universal gas constant (J/mol K)
Patm = 0.101325    # atmospheric pressure (MPa)
### ---------------------------------------------------------------------------


### density (g/cm³) for various fuels
rho_sucrose = 1.58
rho_dextrose = 1.56
rho_fructose = 1.69
rho_sorbitol = 1.49
rho_xylitol = 1.52
rho_erythritol = 1.45
rho_mannitol = 1.52
rho_paraffin = 0.90
rho_charcoal = 1.45

### density (g/cm³) for the oxidizer KNO_3
rho_KNO3 = 2.11

### density for catalyzers
rho_iron_oxide = 5.24
rho_sulphur = 2.0

### densities for structural materials (g/cm³)
rho_steel = 6.777 #7.85
rho_aluminum = 2.65
rho_epoxy = 1.25
rho_PVC = 1.4


#%%% Propellent properties (tables) and choice of fuel type, etc
#%%%% Definition of the Propellant class
class Propellant:
    def __init__(self, prop_name, 
                 oxidizer_name, fuel_name, cata_name, 
                 oxidizer_density, fuel_density, cata_density, 
                 oxi_mass_fraction, fuel_mass_fraction, cata_mass_fraction, 
                 Kn_vs_P_params, r_vs_P_params, pressure_limits, 
                 M, k, To):
        self.prop_name = prop_name       # Name of the propellant, like 'KNSU'
        self.oxidizer_name = oxidizer_name # Name of the oxidizer used in this propellant, like 'potassium nitrate' or 'KNO3'
        self.fuel_name = fuel_name       # Name of the fuel used in this propellant, like 'sucrose'
        self.cata_name = cata_name       # Name of the catalyst used in this propellant, like 'iron oxide' or 'sulphur'
        self.oxidizer_density = oxidizer_density # Density of the oxidizer, in (g/cm³)
        self.fuel_density = fuel_density # Density of the fuel, in (g/cm³)
        self.cata_density = cata_density # Density of the catalyst, in (g/cm³)
        self.oxi_mass_fraction = oxi_mass_fraction # Mass fraction of the oxidizer (ex: 0.65)
        self.fuel_mass_fraction = fuel_mass_fraction # Mass fraction of the fuel (ex: 0.35)
        self.cata_mass_fraction = cata_mass_fraction # Mass fraction of the catalyst (ex: 0)
        self.Kn_vs_P_params = Kn_vs_P_params  # Array for Kv vs P empirical fit parameter values
        self.r_vs_P_params = r_vs_P_params    # Array for r vs P empirical fit parameter values
        self.pressure_limits = pressure_limits  # Array of pressure range limits
        self.M = M  # Effective Molar mass of the products of combustion (g/mol)
        self.k = k  # Specific heats ratio Cp/Cv
        self.To = To  # Ideal combustion temperature (K)

    def display_info(self):
        print(f"Propellant Name: {self.prop_name}")
        print(f"Oxidizer Name: {self.oxidizer_name}")
        print(f"Fuel Name: {self.fuel_name}")
        if self.cata_mass_fraction>0:
            print(f"Catalyst Name: {self.cata_name}")
        print(f"Oxidizer Density: {self.oxidizer_density} g/cm³")
        print(f"Fuel Density: {self.fuel_density} g/cm³")
        if self.cata_mass_fraction>0:
            print(f"Catalyst Density: {self.cata_density} g/cm³")
        print(f"Oxidizer Mass Fraction: {self.oxi_mass_fraction}")
        print(f"Fuel Mass Fraction: {self.fuel_mass_fraction}")
        if self.cata_mass_fraction>0:
            print(f"Catalyst Mass Fraction: {self.cata_mass_fraction}")
        #print(f"Kv vs P Parameters: {self.Kn_vs_P_params}")
        #print(f"r vs P Parameters: {self.r_vs_P_params}")
        #print(f"Pressure Limits: {self.pressure_limits}")
        print(f"Effective Molar Mass (M): {self.M} g/mol")
        print(f"Specific Heats Ratio (k): {self.k}")
        print(f"Ideal Combustion Temperature (To): {self.To} K")

    def calculate_burning_rate(self, P, G, kv):
        for i, limit in enumerate(self.pressure_limits):
            if P < limit:
                r_params = self.r_vs_P_params[i]
                break
        else:
            r_params = self.r_vs_P_params[-1]
        r0, a, n = r_params
        r = (1+kv*G)*a*P**n ### use P in MPa
        return r

    def Kn_vs_P(self, P):
        """ Empirical curve of Kn ('Klemmung') vs pressure (in MPa), for various propellants.
        Polynomial fit up to 6th order, ('params' should have seven entries). """
        #a, b, c, d, e, f, g = Kn_vs_P_params
        Kn = 0.0
        params = self.Kn_vs_P_params
        N_params = len(params)
        for i in range(N_params):
            Kn += params[i]*P**i
        Kn = where(Kn>=0, Kn, 0.0)
        return Kn
    

# Example usage
propellant_example = Propellant(
    prop_name="APCP",
    oxidizer_name="Ammonium Perchlorate",
    fuel_name="Aluminum Powder",
    cata_name="Iron Oxide",
    oxidizer_density=1.95,  # Example density in g/cm³
    fuel_density=2.7,  # Example density in g/cm³
    cata_density=5.2,  # Example density in g/cm³
    oxi_mass_fraction=0.7,  # Example mass fraction
    fuel_mass_fraction=0.2,  # Example mass fraction
    cata_mass_fraction=0.1,  # Example mass fraction
    Kn_vs_P_params=[0.1, 0.2, 0.3],
    r_vs_P_params=[[1e-5, 10.1, 0.6], [1e-5, 7.0, 0.6], [1e-5, 6, 0.5]],
    pressure_limits=[5, 10],  # Example pressure range limits in MPa
    M=26.98,
    k=1.2,
    To=1800  # Example combustion temperature in K
)

# propellant_example.display_info()
# P = 8  # Example pressure in MPa
# G = G_star-(dc**2 -(D0**2-d0**2))/Dt0**2
# burning_rate = propellant_example.calculate_burning_rate(P, G)
# print(f"Burning Rate at {P} MPa: {burning_rate} mm/s")

def get_prop_by_name(name, list_of_props):
    matching_props = []
    for prop in list_of_props:
        if name == prop.prop_name:
            matching_props.append(prop)
    return matching_props

def get_prop_by_fuel(name, list_of_props):
    matching_props = []
    for prop in list_of_props:
        if name == prop.fuel_name:
            matching_props.append(prop)
    return matching_props



        