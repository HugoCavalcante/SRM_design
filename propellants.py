# -*- coding: utf-8 -*-
"""
Objects to be used in my_rocket_motor.py. Mostly proppelants.
Created on Sat Feb 15 13:37:48 2025

@author: hugo
"""


#import rocket_motor_classes 
from rocket_motor_classes import *
list_of_proppelants = []

#%%%% Definitions of the many propellants we can use
knsu = Propellant(
    prop_name="KNSU",
    oxidizer_name="KNO3",
    fuel_name="sucrose",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_sucrose,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[32.954,	44.108,	-1.1025, 0, 0, 0, 0],
    r_vs_P_params=[[1e-5, 8.26, 0.319]],
    pressure_limits=[10.3],  
    M=42.02,
    k=1.133,
    To=1720
)
list_of_proppelants.append(knsu)

kndx = Propellant(
    prop_name="KNDX",
    oxidizer_name="KNO3",
    fuel_name="dextrose",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_dextrose,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[43.500,	0.24168, 50.484, -9.9115, 0, 0, 0],
    r_vs_P_params=[[1e-5, 8.875, 0.619], [0.0, 7.553, -0.009], [0, 3.841, 0.688], [0, 17.20, -0.148], [0, 4.775, 0.442]],
    pressure_limits=[0.779, 2.572, 5.930, 8.502, 11.20],  
    M=42.42,  # effective molar mass of combustion products
    k=1.131,  # ratio of specific heats
    To=1710  # combustion temperature in K
)
list_of_proppelants.append(kndx)

knsb_fine = Propellant(
    prop_name="KNSB fine",
    oxidizer_name="KNO3",
    fuel_name="sorbitol",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_sorbitol,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[34.500,	-10.975, 69.667, -19.664, 2.1478, -8.146E-02, 0],
    r_vs_P_params=[[0, 10.708, 0.625], [0.0, 8.763, -0.314], [0, 7.852, -0.013], [0, 3.907, 0.535], [0, 9.653, 0.064]],
    pressure_limits=[0.807, 1.503, 3.792, 7.033, 10.67],  
    M=39.90,  # effective molar mass of combustion products
    k=1.137,  # ratio of specific heats
    To=1600  # combustion temperature in K
)
list_of_proppelants.append(knsb_fine)

knsb_coarse = Propellant(
    prop_name="KNSB coarse",
    oxidizer_name="KNO3",
    fuel_name="sorbitol",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_sorbitol,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[34.500,	-10.975, 69.667, -19.664, 2.1478, -8.146E-02, 0],
    r_vs_P_params=[[0, 5.130, 0.220]],
    pressure_limits=[10.3],  
    M=39.90,  # effective molar mass of combustion products
    k=1.137,  # ratio of specific heats
    To=1600  # combustion temperature in K
)
list_of_proppelants.append(knsb_coarse)

knxy = Propellant(
    prop_name="KNXY",
    oxidizer_name="KNO3",
    fuel_name="xylitol",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_xylitol,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[34.500,	-10.975, 69.667, -19.664, 2.1478, -8.146E-02, 0],
    r_vs_P_params=[[0, 10.708, 0.625], [0.0, 8.763, -0.314], [0, 7.852, -0.013], [0, 3.907, 0.535], [0, 9.653, 0.064]],
    pressure_limits=[0.807, 1.503, 3.792, 7.033, 10.67],  
    M=39.90,  # effective molar mass of combustion products
    k=1.137,  # ratio of specific heats
    To=1600  # combustion temperature in K
)
list_of_proppelants.append(knxy)

kner = Propellant(
    prop_name="KNER",
    oxidizer_name="KNO3",
    fuel_name="erythritol",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_erythritol,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[50.410, 211.194, -47.7976, 8.779018, -0.830282, 0.030519, 0],
    r_vs_P_params=[[0, 2.9, 0.4]],
    pressure_limits=[10.3],
    M=38.58,  # effective molar mass of combustion products
    k=1.140,  # ratio of specific heats
    To=1608  # combustion temperature in K
)
list_of_proppelants.append(kner)

knfr = Propellant(
    prop_name="KNFR",
    oxidizer_name="KNO3",
    fuel_name="fructose",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_fructose,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[39.748, 51.6180, -0.9264, 0, 0, 0, 0],
    r_vs_P_params=[[0, 7.4, 0.25]],
    pressure_limits=[10.3],
    M=42.42,  # effective molar mass of combustion products
    k=1.131,  # ratio of specific heats
    To=1710  # combustion temperature in K
)
list_of_proppelants.append(knfr)

knmn_coarse = Propellant(
    prop_name="KNMN coarse",
    oxidizer_name="KNO3",
    fuel_name="mannitol",
    cata_name="iron oxide",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_mannitol,  
    cata_density=rho_iron_oxide, 
    oxi_mass_fraction=0.65,  
    fuel_mass_fraction=0.35, 
    cata_mass_fraction=0,  
    Kn_vs_P_params=[38.000,	92.391, -2.0438, 0, 0, 0, 0],
    r_vs_P_params=[[0, 5.130, 0.220]],
    pressure_limits=[10.3],  
    M=39.83,  # effective molar mass of combustion products
    k=1.136,  # ratio of specific heats
    To=1616  # combustion temperature in K
)
list_of_proppelants.append(knmn_coarse)

blackpowder = Propellant(    
    prop_name="Blackpowder",
    oxidizer_name="KNO3",
    fuel_name="charcoal",
    cata_name="sulphur",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_charcoal,  
    cata_density=rho_sulphur, 
    oxi_mass_fraction=0.75,  
    fuel_mass_fraction=0.15, 
    cata_mass_fraction=0.10,  
    Kn_vs_P_params=[50.00, 145.0, 0.0, 0, 0, 0, 0],
    r_vs_P_params=[[0, 20.0, 0.625]],
    pressure_limits=[10.3],  
    M=68.9,  # effective molar mass of combustion products
    k=1.130,  # ratio of specific heats
    To=1770  # combustion temperature in K
)
list_of_proppelants.append(blackpowder)

knpa = Propellant(    
    prop_name="KNPA",
    oxidizer_name="KNO3",
    fuel_name="paraffin",
    cata_name="sulphur",
    oxidizer_density=rho_KNO3,  
    fuel_density=rho_paraffin,
    cata_density=rho_sulphur, 
    oxi_mass_fraction=0.75,  
    fuel_mass_fraction=0.15, 
    cata_mass_fraction=0.10,  
    Kn_vs_P_params=[30.00, 35.0, 0.0, 0, 0, 0, 0],
    r_vs_P_params=[[0, 17.0, 0.625]],
    pressure_limits=[10.3],  
    M=68.9,  # effective molar mass of combustion products
    k=1.130,  # ratio of specific heats
    To=1770  # combustion temperature in K
)
list_of_proppelants.append(knpa)


