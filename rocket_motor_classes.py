# -*- coding: utf-8 -*-
"""
Classes to contain RocketMotor objects

Created on Sat Feb 15 13:30:55 2025

@author: hugo
"""

from numpy import where, pi, radians, tan 
from constants import *

# #%% Constants ----------------------------------------------------------------
# g = 9.81        # standard gravity acceleration
# R_prime = 8.3144598 # universal gas constant (J/mol K)
# Patm = 0.101325    # atmospheric pressure (MPa)
# ### ---------------------------------------------------------------------------


# ### density (g/cm³) for various fuels
# rho_sucrose = 1.58
# rho_dextrose = 1.56
# rho_fructose = 1.69
# rho_sorbitol = 1.49
# rho_xylitol = 1.52
# rho_erythritol = 1.45
# rho_mannitol = 1.52
# rho_paraffin = 0.90
# rho_charcoal = 1.45

# ### density (g/cm³) for the oxidizer KNO_3
# rho_KNO3 = 2.11

# ### density for catalyzers
# rho_iron_oxide = 5.24
# rho_sulphur = 2.0

# ### densities for structural materials (g/cm³)
# rho_steel = 6.777 #7.85
# rho_aluminum = 2.65
# rho_epoxy = 1.25
# rho_PVC = 1.4


class StructuralMaterial:
    def __init__(self, name, density):
        self.name = name
        self.density = density 
    def display_properties(self):
        print('Properties of structural material:')
        print(f'  Name: {self.name}')
        print(f'  Density: {self.density:.3f} g/cm³')
            
        
if __name__ == '__main__':
    ### These are also defined in 'propellants.py'
    steel = StructuralMaterial('steel', rho_steel)
    aluminum = StructuralMaterial('aluminum', rho_aluminum)
    PVC = StructuralMaterial('PVC', rho_PVC)



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

    def display_properties(self):
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


#### We will start a Motor class
### But first, we need a class for each component of the Motor, namelly: Grain, Casing, Bulkhead and Nozzle

class Grain:
    def __init__(self, Dg0, dg0, Lg0, Nsegs, Lseg, inhibit_ext, propellant, hc, G_star, kv):
        self.Dg0 = Dg0  # External diameter of the grain (initial)
        self.dg0 = dg0  # Diameter of the grain core (initial)
        self.Lg0 = Lg0  # Initial grain segment length (total)
        self.Nsegs = Nsegs  # Number of grain segments
        self.Lseg = Lseg  # Initial length of each segment
        self.inhibit_ext = inhibit_ext  # Boolean to indicate if external surface is inhibited from burning
        self.propellant = propellant  # Propellant used in the grain
        self.hc = hc  # Combustion efficiency
        self.G_star = G_star  # Propellant 'erosive burning' area ratio threshold
        self.kv = kv  # Propellant 'erosive burning' velocity coefficient

    def grain_volume(self):
        """ grain volume, in cm³ """
        return pi*((self.Dg0/10)**2-(self.dg0/10)**2)*self.Lg0/40

    def grain_mass(self):
        rho_p = 1.0/(
            self.propellant.oxi_mass_fraction/self.propellant.oxidizer_density + 
            self.propellant.fuel_mass_fraction/self.propellant.fuel_density +
            self.propellant.cata_mass_fraction/self.propellant.cata_density)
        non_ideal_density_factor = 0.95
        rho_p = rho_p*non_ideal_density_factor
        return self.grain_volume()*rho_p

    def display_properties(self):
        print('Grain Properties:')
        print(f'  External Diameter (Dg0): {self.Dg0} mm')
        print(f'  Core Diameter (dg0): {self.dg0} mm')
        print(f'  Total Length (Lg0): {self.Lg0:.1f} mm')
        print(f'  Number of Segments (Nsegs): {self.Nsegs}')
        print(f'  Length of Each Segment (Lseg): {self.Lseg:.1f} mm')
        print(f'  Inhibit External Surface: {self.inhibit_ext}')
        print(f'  Initial Propellant Mass: {self.grain_mass():.0f} g')
        print(f'  Combustion Efficiency (hc): {self.hc}')
        print(f'  Erosive Burning Area Ratio Threshold (G_star): {self.G_star}')
        print(f'  Erosive Burning Velocity Coefficient (kv): {self.kv}')
        print(f'  Propellant name: {self.propellant.prop_name}')

# Example of creating a Grain object
if __name__ == "__main__":
    grain = Grain(Dg0=46, dg0=30, Lg0=300, Nsegs=1, Lseg=(300-2)/1, inhibit_ext=True, propellant = propellant_example, hc=0.95, G_star=2.5, kv=0.05)
    grain.display_properties()


class CombustionChamber:
    def __init__(self, Dc, dc, Lc, Bo, No, material):
        self.Dc = Dc  # Chamber external diameter (in mm)
        self.dc = dc  # Chamber internal diameter (in mm)
        self.Lc = Lc  # Chamber length (in mm)
        self.Bo = Bo  # Bulkhead offset (in mm)
        self.No = No  # Nozzle offset (in mm)
        self.t = (Dc-dc)/2
        self.material = material  # Material of the combustion chamber

    def chamber_volume(self):
        """ chamber available (empty) volume, in mm³ """
        return pi*(self.dc**2)*self.Lc/4

    def casing_volume(self):
        """ casing volume, in cm³ """
        return pi*((self.Dc/10)**2-(self.dc/10)**2)*(self.Lc+self.Bo+self.No)/40

    def casing_mass(self):
        return self.casing_volume()*self.material.density

    def display_properties(self):
        print('Combustion Chamber Properties:')
        print(f'  External Diameter (Dc): {self.Dc} mm')
        print(f'  Internal Diameter (dc): {self.dc} mm')
        print(f'  Chamber Length (Lc): {self.Lc} mm')
        print(f'  Bulkhead Offset (Bo): {self.Bo} mm')
        print(f'  Nozzle Offset (No): {self.No} mm')
        print(f'  Chamber Volume (Empty): {self.chamber_volume():.0f} mm³')
        print(f'  Casing Mass: {self.casing_mass():.0f} g')
        print(f'  Material: {self.material.name}')

# Example of creating a CombustionChamber object
if __name__ == "__main__":
    chamber = CombustionChamber(Dc=49, dc=46, Lc=302, Bo=30, No=30, material=steel)
    chamber.display_properties()


class Bulkhead:
    def __init__(self, material, thickness, mass):
        self.material = material
        self.thickness = thickness
        self.mass = mass
    
    def display_properties(self):
        print('Bulkhead Properties:')
        print(f'  Material: {self.material.name}')
        print(f'  Thickness: {self.thickness} mm')
        print(f'  Mass: {self.mass} g')
    
if __name__ == "__main__":
    bulkhead = Bulkhead(material = steel, thickness=25, mass=50)

class Nozzle:
    def __init__(self, D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0):
        self.D1 = D1  # Entry diameter (in mm)
        self.Dt0 = Dt0  # Initial throat diameter (in mm)
        self.Dtf = Dtf if Dtf is not None else Dt0  # Final throat diameter (in mm, defaults to Dt0)
        self.Lt = Lt  # Length of the throat (in mm, defaults to 0)
        self.beta = beta  # Convergence half-angle (in degrees)
        self.alpha = alpha  # Divergence half-angle (in degrees)
        self.De = De if De is not None else D1  # Exit diameter (in mm, defaults to D1)
        self.mass = mass # Nozzle mass (g)

    def calculate_convergent_length(self):
        beta_rad = radians(self.beta)
        L_nozzle_c = ((self.D1 - self.Dt0) / 2) / tan(beta_rad)
        return L_nozzle_c

    def calculate_divergent_length(self):
        alpha_rad = radians(self.alpha)
        L_nozzle_d = ((self.De - self.Dt0) / 2) / tan(alpha_rad)
        return L_nozzle_d

    def calculate_total_length(self):
        L_nozzle_c = self.calculate_convergent_length()
        L_nozzle_d = self.calculate_divergent_length()
        L_nozzle = L_nozzle_c + L_nozzle_d + self.Lt
        return L_nozzle


    def draw_nozzle(self):
        if not "plot" in globals():
            from matplotlib.pyplot import plot
        if not "subplots" in globals():
            from matplotlib.pyplot import subplots 
        if not "grid" in globals():
            from matplotlib.pyplot import grid
        if not "arrow" in globals():
            from matplotlib.pyplot import arrow 
        L_nozzle_c = self.calculate_convergent_length()
        L_nozzle_d = self.calculate_divergent_length()
        L_nozzle = self.calculate_total_length()
        alpha_rad = radians(self.alpha)
        beta_rad = radians(self.beta)
        fig, ax = subplots()
        
        plot([-self.Lt/2-(self.Dtf-self.Dt0)/tan(alpha_rad)/2, self.Lt/2+(self.Dtf-self.Dt0)/tan(beta_rad)/2],
             [self.Dtf/2, self.Dtf/2], lw = 2, color = 'red')
        plot([-self.Lt/2-(self.Dtf-self.Dt0)/tan(alpha_rad)/2, self.Lt/2+(self.Dtf-self.Dt0)/tan(beta_rad)/2],
             [-self.Dtf/2, -self.Dtf/2], lw = 2, color = 'red')
        plot([-L_nozzle_d-self.Lt/2, -self.Lt/2, self.Lt/2, self.Lt/2+L_nozzle_c], 
             [self.De/2, self.Dt0/2, self.Dt0/2, self.D1/2], lw = 2, color = 'black')
        plot([-L_nozzle_d-self.Lt/2, -self.Lt/2, self.Lt/2, self.Lt/2+L_nozzle_c], 
             [-self.De/2, -self.Dt0/2, -self.Dt0/2, -self.D1/2], lw = 2, color = 'black')
        
        ax_x_lims = ax.get_xlim()
        ax_y_lims = ax.get_ylim()
        ax.set_box_aspect((ax_y_lims[1]-ax_y_lims[0])/(ax_x_lims[1]-ax_x_lims[0]))
        grid()
        arrow(0,0, -L_nozzle/5, 0, width = 0.8)


    def display_properties(self, draw_nozzle=True):
        print('Nozzle Properties:')
        print(f'  Entry Diameter (D1): {self.D1} mm')
        print(f'  Initial Throat Diameter (Dt0): {self.Dt0} mm')
        print(f'  Final Throat Diameter (Dtf): {self.Dtf} mm')
        print(f'  Length of the Throat (Lt): {self.Lt} mm')
        print(f'  Convergence Half-Angle (beta): {self.beta} degrees')
        print(f'  Divergence Half-Angle (alpha): {self.alpha} degrees')
        print(f'  Exit Diameter (De): {self.De} mm')
        print(f'  Convergent Length: {self.calculate_convergent_length():.2f} mm')
        print(f'  Divergent Length: {self.calculate_divergent_length():.2f} mm')
        print(f'  Total Length: {self.calculate_total_length():.2f} mm')
        print(f'  Mass: {self.mass} g')
        if draw_nozzle:
            self.draw_nozzle()

# Example of creating a Nozzle object
if __name__ == "__main__":
    nozzle = Nozzle(D1=46, Dt0=15, beta=30, alpha=12, mass = 70)
    nozzle.display_properties()


### Finally, the class for the motor
class RocketMotor:
    def __init__(self, name, grain, combustion_chamber, bulkhead, nozzle):
        self.name = name 
        self.grain = grain  # Grain object
        self.combustion_chamber = combustion_chamber  # CombustionChamber object
        self.bulkhead = bulkhead  # Bulkhead object
        self.nozzle = nozzle  # Nozzle object
        self.empty_mass = self.bulkhead.mass+self.combustion_chamber.casing_mass()+self.nozzle.mass
        self.initial_mass = self.empty_mass + self.grain.grain_mass()

    def display_properties(self):
        print('Rocket Motor Properties:')
        print(f'  Motor name: {self.name}')
        #self.grain.display_properties()
        #self.combustion_chamber.display_properties()
        #self.bulkhead.display_properties()
        #self.nozzle.display_properties()
        motor_length = (self.combustion_chamber.Lc + self.combustion_chamber.Bo + self.nozzle.calculate_total_length())
        #volumetric_loading_fraction
        VLF = 1000*self.grain.grain_volume()/self.combustion_chamber.chamber_volume()
        #Ap = (pi/4)*self.grain.dg0**2
        #At = (pi/4)*self.nozzle.Dt0**2
        #PtTR = Ap/At
        PtTR = (self.grain.dg0/self.nozzle.Dt0)**2
        print(f'  Empty Mass: {self.empty_mass:.0f} g')
        print(f'  Propellant Mass (mg): {self.grain.grain_mass():.0f} g')
        print(f'  Initial Mass: {self.initial_mass:.0f} g')
        print(f'  Diameter (Dc): {self.combustion_chamber.Dc:.0f} mm')
        print(f'  Chamber Internal Length: {self.combustion_chamber.Lc:.0f} mm')
        print(f'  Length: {motor_length:.0f} mm')
        print(f'  Throat Diameter (Dt): {self.nozzle.Dt0:.0f} mm')
        print(f'  Volumetric Loading Fraction (VLF): {100*VLF:.1f}%')
        print(f'  Port-to-throat ratio (Ap/At): {PtTR:.1f}')


# Example of creating a RocketMotor object
if __name__ == "__main__":
    rocket_motor = RocketMotor("Example Motor", grain, chamber, bulkhead, nozzle)
    rocket_motor.display_properties()


def get_motor_by_name(name, list_of_motors):
    matching_motors = []
    for motor in list_of_motors:
        if name == motor.name:
            matching_motors.append(motor)
    return matching_motors



        