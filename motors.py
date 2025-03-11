# -*- coding: utf-8 -*-
"""

Objects to be used in my_rocket_motor.py. Motors and their components (grains, casing, bulkhead, nozzle).

Created on Mon Mar 10 10:53:33 2025
@author: hugo
"""

from constants import *
from rocket_motor_classes import *
from propellants import *

list_of_motors = []
list_of_grains = []
list_of_casings = []
list_of_bulkheads = []
list_of_nozzles = []

###############################################################################
#### Richard Nakka's Impulser (2013), class I, four BATES segments with inhibited outer surface.
#### May use KNDX or KNSB
D0 = 34.8 # mm (1.5" - 2*0.065")
d0 = 9.525 # mm (3/8")
Nsegs = 4
L0 = Nsegs*(1/2)*(3*D0+d0) # 227.85 mm
Lseg = L0/Nsegs
Dc = 38.1 # mm (1.5")
dc = D0
Lc = 293.116 # mm (11.54")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 7.053 #6.706 # mm (0.264")
Dtf = 7.086 #6.731
Lt = 1 # mm (length of the throat region)
inhibit_ext = True
beta = 30 # nozzle covergence angle (degrees)
alpha = 10 # nozzle divergence angle (degrees)
De = dc # exit diameter at the end of the nozzle (mm)
## creating the objects
impulser_knsb_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = knsb, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(impulser_knsb_grain)

impulser_bulkhead = Bulkhead(aluminum, 20, rho_aluminum*2*pi*(d0**2)/4)
list_of_bulkheads.append(impulser_bulkhead)

impulser_casing = CombustionChamber(Dc, dc, Lc, 20, 20, aluminum)
list_of_casings.append(impulser_casing)


### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
impulser_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(impulser_nozzle)

impulser_knsb = RocketMotor("Impulser KNSB", grain = impulser_knsb_grain, combustion_chamber = impulser_casing, bulkhead=impulser_bulkhead, nozzle=impulser_nozzle)
list_of_motors.append(impulser_knsb)
####

#######################################

#### Impulser with KNDX
impulser_kndx_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = kndx, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(impulser_kndx_grain)

impulser_kndx_nozzle = Nozzle(dc, 7.728, 7.761, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(impulser_kndx_nozzle)

impulser_kndx = RocketMotor("Impulser KNDX", grain = impulser_kndx_grain, combustion_chamber = impulser_casing, bulkhead=impulser_bulkhead, nozzle=impulser_kndx_nozzle)
list_of_motors.append(impulser_kndx)
####
###############################################################################

###############################################################################
#### Impulser X (longer casing, five segments)
D0 = 34.8 # mm (1.5" - 2*0.065")
d0 = 9.525 # mm (3/8")
Nsegs = 5
L0 = Nsegs*(1/2)*(3*D0+d0) # 284.8 mm
Lseg = L0/Nsegs
Dc = 38.1 # mm (1.5")
dc = D0
Lc = 356.6 # mm (14.04")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 8.640 #8.763 # mm (0.345")
Dtf = 8.720 #8.843 
Lt = 1 # mm (length of the throat region)
inhibit_ext = True
beta = 30 # nozzle covergence angle (degrees)
alpha = 10 # nozzle divergence angle (degrees)
De = dc # exit diameter at the end of the nozzle (mm)

## creating the objects
impulserX_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = kndx, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(impulserX_grain)

impulserX_casing = CombustionChamber(Dc, dc, Lc, 20, 20, aluminum)
list_of_casings.append(impulserX_casing)
### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
impulserX_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(impulserX_nozzle)
impulserX = RocketMotor("Impulser X", grain = impulserX_grain, combustion_chamber = impulserX_casing, bulkhead=impulser_bulkhead, nozzle=impulserX_nozzle)
list_of_motors.append(impulserX)
#########
###############################################################################

###############################################################################
#### Impulser XX (even longer casing, six segments)
D0 = 34.8 # mm (1.5" - 2*0.065")
d0 = 9.525 # mm (3/8")
Nsegs = 6
L0 = Nsegs*(1/2)*(3*D0+d0) # 
Lseg = L0/Nsegs
Dc = 38.1 # mm (1.5")
dc = D0
Lc = 420.37 # mm (16.55")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 9.804 #8.763 # mm (0.345")
Dtf = 9.884 #8.843 
inhibit_ext = True
Lt = 1 
beta = 30 # nozzle covergence angle (degrees)
alpha = 10 # nozzle divergence angle (degrees)
De = dc # exit diameter at the end of the nozzle (mm)
## creating the objects
impulserXX_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = kndx, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(impulserXX_grain)

impulserXX_casing = CombustionChamber(Dc, dc, Lc, 20, 20, aluminum)
list_of_casings.append(impulserXX_casing)

### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
impulserXX_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(impulserXX_nozzle)

impulserXX = RocketMotor("Impulser XX", grain = impulserXX_grain, combustion_chamber = impulserXX_casing, bulkhead=impulser_bulkhead, nozzle=impulserXX_nozzle)
list_of_motors.append(impulserXX)


###
##### Design parameter for the Impulser et al., according to Nakka
###___________________________________________________________________________
### Parameter      | units |   Impulser        |  Impulser-X  |  Impulser-XX  |
### Propellant     |   -   |   KNSB   |  KNDX  |    KNDX      |    KNDX       |
### Propel. mass   |   g   |   300    |   306  |    383       |     465       |
### MEOP (max P)   |  MPa  |   6.9    |   6.9  |    6.9       |     6.9       |
### Kn range       |   -   | 324-410  | 270-342|   270-342    |   265-330     |
### throat diam.   |  mm   |   7.053  | 7.728  |    8.640     |    9.804      |
### Max thrust     |   N   |    400   |   500  |     625      |     750       |
### Thrust duration|   s   |    1.1   |   0.9  |     0.9      |     0.9       |
### Total impulse  |   Ns  |    380   |   400  |     500      |     600       |
###________________|_______|__________|________|______________|_______________|
###
###
###############################################################################

###############################################################################
#### Richard Nakka's Kappa (2001), class K, four BATES segments with inhibited outer surface.
D0 = 55.1 # mm (2.170")
d0 = 19.05# mm (3/4")
Nsegs = 4
L0 = Nsegs*(1/2)*(3*D0+d0) # 
Lseg = L0/Nsegs
Dc = 63.5 # mm (2.5")
dc = D0
Lc = 462 # mm (18.2")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 12.75 # mm (0.502")
Dtf = 12.75 
Lt = 1 
inhibit_ext = True
beta = 45 # nozzle covergence angle (degrees)
alpha = 12 # nozzle divergence angle (degrees)
De = 53.85 # exit diameter at the nozzle (mm)
propellant = kndx 

## creating the objects
kappa_kndx_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = propellant, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(kappa_kndx_grain)

kappa_bulkhead = Bulkhead(steel, 30, rho_steel*3*pi*(d0**2)/4)
list_of_bulkheads.append(kappa_bulkhead)

kappa_casing = CombustionChamber(Dc, dc, Lc, 30, 30, steel)
list_of_casings.append(kappa_casing)

### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
kappa_kndx_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(kappa_kndx_nozzle)

kappa_kndx = RocketMotor("Kappa KNDX", grain = kappa_kndx_grain, combustion_chamber = kappa_casing, bulkhead=kappa_bulkhead, nozzle=kappa_kndx_nozzle)
list_of_motors.append(kappa_kndx)
### It would be nice to have the __init__() function in those classes to auto append a created object in a global list of existing objects. I guess it can be made later, with dir() and for obj in dir if type(obj)==something ...
####

### Results for the Kappa
###____________________________________________
### Parameter      | units |       Kappa       |
### Propellant     |   -   |   KNDX   |  KNSB  |
### Propel. mass   |   g   |   1500   |  1500  |
### MEOP (max P)   |  MPa  |   8.5    |   8.5  |
### Kn range       |   -   | 320-380  |    -   |
### throat diam.   |  mm   |   12.75  |    -   |
### Max thrust     |   N   |    1580  |  1620  |
### Thrust duration|   s   |    1.5   |1.5-2.0 |
### Total impulse  |   Ns  |   2003   |  1821* |
### Specif. impulse|   s   |   137    |  125   |
###________________|_______|__________|________|
#### *predicted 1987 Ns
###
###############################################################################

###############################################################################
#### Richard Nakka's Juno, class J, single hollow grain, no inhibitors, regressive profile
D0 = 45 # mm (2.170")
d0 = 15 # ?  
Nsegs = 1
#L0 = Nsegs*(1/2)*(3*D0+d0) # 
L0 = 300 # ?
Lseg = L0/Nsegs
Dc = 48 # mm (1.875")
dc = D0
Lc = 338 # mm (13.3")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 15 # mm (0.590")
Dtf = 15 
inhibit_ext = False
Lt = 1 # mm (throat length)
beta = 30 # nozzle covergence angle (degrees)
alpha = 12 # nozzle divergence angle (degrees)
De = 42.4 # exit diameter at the nozzle (mm)
propellant = kndx 

## creating the objects
juno_kndx_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = propellant, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(juno_kndx_grain)

juno_bulkhead = Bulkhead(steel, 30, rho_steel*3*pi*(d0**2)/4)
list_of_bulkheads.append(juno_bulkhead)

juno_casing = CombustionChamber(Dc, dc, Lc, 30, 30, steel)
list_of_casings.append(juno_casing)

### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
juno_kndx_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=dc, mass = 30)
list_of_nozzles.append(juno_kndx_nozzle)

juno_kndx = RocketMotor("Juno KNDX", grain = juno_kndx_grain, combustion_chamber = juno_casing, bulkhead=juno_bulkhead, nozzle=juno_kndx_nozzle)
list_of_motors.append(juno_kndx)


### Results for the Juno
###____________________________________________
### Parameter      | units |       Juno        |
### Propellant     |   -   |   KNDX   |  KNSB  |
### Propel. mass   |   g   |    650   |        |
### MEOP (max P)   |  MPa  |    6.9   |        |
### Kn range       |   -   |  285-312 |    -   |
### throat diam.   |  mm   |    15    |    -   |
### Max thrust     |   N   |    1740  |        |
### Thrust duration|   s   |    0.56  |        |
### Total impulse  |   Ns  |    885   |        |
### Specif. impulse|   s   |    139   |        |
###________________|_______|__________|________|
###
######
###############################################################################

###############################################################################
##### OBAFOG motor (PVC), no nozzle
D0 = 17 # mm
d0 = 5 # mm 
Nsegs = 1
#L0 = Nsegs*(1/2)*(3*D0+d0) # 
L0 = 94 
Lseg = L0/Nsegs
Dc = 20 # mm 
dc = 17
Lc = 100 # mm (13.3")
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 7 # mm 
Dtf = 16 
inhibit_ext = False
Lt = 1.5 # mm (throat length)
beta = 85 # nozzle covergence angle (degrees)
alpha = 87 # nozzle divergence angle (degrees)
De = 7.1 # exit diameter at the nozzle (mm)
propellant = knsu 

## creating the objects
OBA_grain = Grain(Dg0=D0, dg0=d0, Lg0=L0, Nsegs=Nsegs, Lseg=Lseg, inhibit_ext=inhibit_ext, propellant = propellant, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(OBA_grain)

OBA_bulkhead = Bulkhead(PVC, 3, rho_PVC*0.3*pi*(d0**2)/4)
list_of_bulkheads.append(OBA_bulkhead)

OBA_casing = CombustionChamber(Dc, dc, Lc, 2, 2, PVC)
list_of_casings.append(OBA_casing)

### Nozzle: Nozzle(D1, Dt0, Dtf=None, Lt=0, beta=0, alpha=0, De=None, mass=0)
OBA_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De=De)
list_of_nozzles.append(OBA_nozzle)

OBA_motor = RocketMotor("OBAFOG KNSU", grain = OBA_grain, combustion_chamber = OBA_casing, bulkhead=OBA_bulkhead, nozzle=OBA_nozzle)
list_of_motors.append(OBA_motor)
###############################################################################

###############################################################################
##### Thais 1 - Variant of OBAFOG model
Thais_grain = Grain(Dg0= 25.4-2, dg0=7, Lg0=42, Nsegs=1, Lseg=42, inhibit_ext=False, propellant=knxy, hc=0.95, G_star=1, kv=0.03)
list_of_grains.append(Thais_grain)

Thais_bulkhead = Bulkhead(PVC, 3, rho_PVC*0.3*pi*(d0**2)/4)
list_of_bulkheads.append(Thais_bulkhead)

Thais_casing = CombustionChamber(Dc=25.4, dc=25.4-2.0, Lc=47.0, Bo=1, No=1, material=PVC)
list_of_casings.append(Thais_casing)

Thais_nozzle = Nozzle(D1=23.4, Dt0=5, Dtf=14, beta=87, alpha=87, De=5, mass=12)
list_of_nozzles.append(Thais_nozzle)

Thais_1 = RocketMotor(name = "Thais 1", grain=Thais_grain, combustion_chamber=Thais_casing, bulkhead=Thais_bulkhead, nozzle=Thais_nozzle)
list_of_motors.append(Thais_1)
###############################################################################


###############################################################################
##### Carcara - v0 (steel motor, originally thought to reach 1 km)
D0 = 42 # mm
d0 = 26 # mm 
Nsegs = 1
#L0 = Nsegs*(1/2)*(3*D0+d0) # 
L0 = 294 
Lseg = L0/Nsegs
Dc = 48.5 # mm 
dc = 44.5
Lc = 300 # mm 
Bo = 14 # mm (offset for the bulkhead: L_ext = Lc + B0 + No)
No = 14 # mm (offset for the nozzle)
thickness = (Dc - dc)/2  # thickness of the combustion chamber wall (mm)
Dt0 = 18 # mm 
Dtf = 18 
inhibit_ext = False
Lt = 2.0 # mm (throat length)
beta = 85 # nozzle covergence angle (degrees)
alpha = 85 # nozzle divergence angle (degrees)
De = 18.1 # exit diameter at the nozzle (mm)
propellant = knxy 

carcara_grain = Grain(D0, d0, L0, Nsegs, Lseg, inhibit_ext, propellant, 0.95, 1.0, 0.03)
list_of_grains.append(carcara_grain)
carcara_bulkhead = Bulkhead(steel, 18, 108)
list_of_bulkheads.append(carcara_bulkhead)
carcara_casing = CombustionChamber(Dc, dc, Lc, Bo, No, steel)
### m_casing is 594 g (measured)
### total motor mass (empty) should be 752 g
list_of_casings.append(carcara_casing)
carcara_nozzle = Nozzle(dc, Dt0, Dtf, Lt, beta, alpha, De, 50)
carcara = RocketMotor("Carcara v0", carcara_grain, carcara_casing, carcara_bulkhead, carcara_nozzle)
list_of_motors.append(carcara)
###############################################################################




