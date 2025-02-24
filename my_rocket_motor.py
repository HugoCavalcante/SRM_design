# -*- coding: utf-8 -*-
#%% File header - my_rocket_motor.py
"""
Calculates the thrust of experimental rocket motors.
Based on Richard Nakka's Solid Rocket Motor Theory at https://www.nakka-rocketry.net/th_intro.html
Created on Thu Feb  6 11:19:46 2025

@author: hugo
"""

#%% Modules imports
from numpy import linspace, zeros, zeros_like, pi, sqrt, log2, mean, where, ceil, tan
from matplotlib.pyplot import figure, legend, show, xlabel, ylabel, plot
from matplotlib.pyplot import subplot, subplots, tick_params, title
from matplotlib.pyplot import tight_layout, grid, arrow
from rocket_motor_classes import *
from propellants import *

#%% Choose which graphs to plot
plot_Kn_vs_x = True 
plot_Kn_vs_P = True
plot_grain_geometry = True
plot_P_vs_t = True
plot_combustion_vs_t = True
plot_F_vs_t = True
draw_nozzle = True
plot_nozzle_vars = True

### A few parameters affecting the burning rate
hc = 0.95       # combustion efficiency
G_star = 1.0    # propellant 'erosive burning' area ratio threshold
kv = 0.03       # propellant 'erosive burning' velocity coefficient


#%%%  Parameters -------------------------------------------------------------
######### Propellant grain parameters
# Grain external diameter
D0 = 46.0 #21.4 #44.0 #369.0 #43.0 # External diameter of propellant grain, in mm (initial)
d0 = 30.0 #29.0 #18.0 # Diameter of the propellant grains hollow core, in mm (initial)
L0 = 294 #48.0 #294.0 #290.0 # Length of propellant grain (full, all segments), in mm
Nsegs = 1 #2 #4 #1 # Number of segments comprising the propellant grain
Lseg = L0/Nsegs # Length of each segment of the propellant grain

### Inhibit external cylinder surface (for the BATES method)
inhibit_ext = False

######### Combustion chamber parameters
Dc = 48.5 #25.4 #48.5 #75.0 #48.5  # Internal diameter of the chamber (mm)
Lc = 300 #48.0 #300.0 #470.0 #290.0 # Internal chamber length  (mm)
### You can specify eiher dc or thickness
###---------------------------------------------
#dc = 44.5 # Internal chamber diameter (mm)
#thickness = (Dc - dc)/2  # thichness of the combustion chamber wall (mm)
thickness = 2.0 #1.8 # thichness of the combustion chamber wall (mm)
dc = Dc - 2*thickness
###---------------------------------------------
### throat diameter (mm) (initial)
Dt0 = 18.0 #5.0 #15.6 #10.0 #18.0
### throat diameter (mm) (final) 
Dtf = 20.0 #18.0 #14.0
### Rate of increase in throat diameter (due to erosion) to regression. Usually set to 0, but important for a hole in a PVC cap. 
e_rate = 2*(Dtf-Dt0)/(D0-d0)

mass_casing_wall = 594 #3.0 #594.0 # grams (just the tube) 
mass_bulkhead = 108 #0.0 #108.0 #52.0 # grams
mass_nozzle = 50 #5.0 #50.0 #40.0   # grams
mc = mass_casing_wall + mass_bulkhead + mass_nozzle 


#%%%% Choice to the propellant to be used
propellant = knsu
### You can see the available propellants in the variable 'list_of_propellants'
# print(list_of_propellants)
### You can also get the propellant object by its exact name:
# propellant = get_prop_by_name("KNSB fine", list_of_propellants)
### or a list of propellants made with a specific fuel
#print(get_prop_by_name("sorbitol", list_of_propellants)) 

prop_name, fuel_name = propellant.prop_name, propellant.fuel_name
rho_f = propellant.fuel_density
Kn_vs_P_params = propellant.Kn_vs_P_params

print(f"Using propellant {prop_name} (KNO3 + {fuel_name})")
propellant.display_info()

### Mass densities of oxidizer rho_o and fuel rho_f (g/cm³)
# rho_o = rho_KNO3
# rho_f = rho_sucrose
# rho_cata = rho_sulphur
rho_o = propellant.oxidizer_density
rho_f = propellant.fuel_density
rho_cata = propellant.cata_density


### Propellent composition (mass fractions of oxidizer, fuel and catalyst)
### A standard is 65/35/0
### For blackpowder: 75% oxi/15% charcoal / 10% sulphur
mass_fraction_oxi = propellant.oxi_mass_fraction
mass_fraction_fuel = propellant.fuel_mass_fraction
mass_fraction_cata = propellant.cata_mass_fraction

### propellant mass density (g/cm³)
rho_p = 1.0/(mass_fraction_oxi/rho_o + mass_fraction_fuel/rho_f +mass_fraction_cata/rho_cata)
non_ideal_density_factor = 0.95
rho_p = rho_p*non_ideal_density_factor

print(f"Propellent effective density = {rho_p:.4g} g/cm³.")

### Choosing the material for the combustion chamber
rho_chamber = rho_steel


#%% Functions definitions ----------------------------------------------------
# def flowrate(Ab, r):
#     """ Mass flow-rate equation. Right-hand side of the differential equation for dm/dt.
#         dm_g/dt = A_b rho_p r, where m_g is the ejected mass, A_b is the burning area,
#         and r is the burning rate (length/time)."""
#     return Ab*rho_p*r

def hcylinder_area(D, d, L):
    """ Area of the surfaces of a hollow cylinder.
    A = pi*D*L + pi*d*L + pi(D²-d²)/2. D is the external diameter, d is the hollow diameter,
    L is the length (height). """
    return pi*(L*(D+d) +(D**2-d**2)/2)

def BATES_area(D, d, L):
    """ Area of the surfaces of a hollow cylinder  without external side surface (BATES).
    A = pi*d*L + pi(D²-d²)/2. D is the external diameter, d is the hollow diameter,
    L is the length (height). """
    return pi*(L*d +(D**2-d**2)/2)

def cyl_vol(D, d, L):
    """ Volume of a hollow cylinder."""
    return pi*(D**2-d**2)*L/4

def deSaintRobert(P):
    """ Model of de Saint Robert for burn rate dependence on pressure.
        r = r_0 a*Pc^n , r_0, a and n are constants. """
    r_0 = 0; a = 7.0; n = 0.625
    return r_0 +a*P**n

# def Kn_vs_P(P, params):
#     """ Empirical curve of Kn ('Klemmung') vs pressure (in MPa), for various propellants.
#     Polynomial fit up to 6th order, ('params' should have seven entries). """
#     #a, b, c, d, e, f, g = params
#     Kn = 0.0
#     N_params = len(params)
#     for i in range(N_params):
#         Kn += params[i]*P**i
#     return Kn

### ---------------------------------------------------------------------------


#%% Initial calculations
### Check for errors in the specified parameters
if (Dc < D0) or (Lc < L0):
    print("Your grain will not fit the casing. Check the parameters!")
    raise SystemExit
elif D0 <= d0:
    print("Your grain core is too big! Check the parameters!")
    raise SystemExit
elif Dc <= dc:
    print("The internal diameter of the chamber dc should be smaller than the external diameter Dc! Check the parameters!")
    raise SystemExit


### Initial mass of grain in grams, (estimated from geometry and density)
# Volume of propellant (cm³)
Vp = cyl_vol(D0/10, d0/10, L0/10) # Notice that we transform the dimensions given in mm to cm
mg0 = rho_p*Vp

### we can estimate the mass of the casing mc, if not known
##mass_casing_wall = rho_chamber*cyl_vol(Dc/10, dc/10, Lc/10)
mc = mass_casing_wall + mass_bulkhead + mass_nozzle
m0 = mg0 + mc

print(f"Mass prop = {mg0:.4g} g, Mass casing = {mc:.4g} g, Total mass = {m0:.4g} g.")

Va = cyl_vol(dc/10, 0.0, Lc/10) # volume available in the chamber (cm³)
### Volumetric loading fraction (Vl) 
Vl = Vp/Va
### "Web thickness" Wf = (Dp - dp)/Dp
Wf0 = (D0-d0)/D0
### Burning area (mm²)
Ab0 = Nsegs*hcylinder_area(D0, d0, L0/Nsegs)
### Throat area (mm²)
At0 = pi*(Dt0**2)/4
Ap0 = pi*(D0**2)*(1-Vl)/4

print(f"Volumetric loading fraction = {Vl:.4g}, Web fraction  = {Wf0:.4g}, Port-to-throat ratio = {Ap0/At0:.4g}.")

#%% Kn vs P
### Kn ('Klemmung') is the ratio Ab/At and is modeled as depending on the chamber pressure.

#Kn_vs_P_params = set_Kn_vs_P_params(Fuel_name)
#Kn_vs_P_params = propellant.Kn_vs_P_params

### Block to plot Kn vs P for the given fuel type
#print(f"Using propellant {fuel_name}")
if plot_Kn_vs_P:
    figure('Kn_vs_P_model')
    P_plot = linspace(0, 10, num = 3000)
    #plot(P_plot, Kn_vs_P(P_plot, Kn_vs_P_params), label = prop_name)
    plot(P_plot, propellant.Kn_vs_P(P_plot), label = prop_name)
    xlabel("$P$ (MPa)")
    ylabel("$Kn$")
    legend()
###

#%% Kn vs x (web regression)
### Building the first table in Nakka's spreadsheet (SRM_2023) Kn vs web regression
N_web = 30000  # number of points
x0 = 0.0      # initial regression
if inhibit_ext:
    xf = min( (D0-d0)/2 , L0/2/Nsegs)   # final regression
else:
    xf = min ( (D0-d0)/4 , L0/2/Nsegs)   # final regression
#xf = (D0-d0)/2  # final regression
dx = (xf -x0)/(N_web-1)  # regression step
x = linspace(x0, xf, num = N_web)  # regression of the burning surfaces (mm)
d = d0+2*x    # internal diameter
if inhibit_ext:
    D = D0   # external diameter (constant if inhibited)
else:
    D = D0 -2*x # external diameter (regressing if not inhibited)
L = L0 - 2*Nsegs*x  # total length regresses in every exposed face
tweb = (D-d)/2      # Web thickness (nonnormalized, mm)

if inhibit_ext:
    Ab = Nsegs*BATES_area(D, d, L/Nsegs)
else:
    Ab = Nsegs*hcylinder_area(D, d, L/Nsegs)

#Dt = Dt0+e_rate*(tweb[0]-tweb)/tweb[0]
Dt = Dt0+e_rate*(tweb[0]-tweb)
At = (pi/4)*Dt**2
Kn = Ab/At
print(f"Kn max = {max(Kn):.4g}, Kn min = {min(Kn):.4g}, Kn avg = {mean(Kn):.4g}")

#%%%
### Block to plot Kn vs web regression (x) ------------------------------------
if plot_Kn_vs_x:
    #print(f"Using propellant {fuel_name}")
    figure('Kn_vs_x')
    plot(x, Kn, label = prop_name)
    ylabel("Kn")
    xlabel("Web regression $x$ (mm)")
    legend()

if plot_grain_geometry:
    figure('Grain geometry')
    if not(inhibit_ext):
        plot(x, D, label = 'D')
    else:
        plot([x[0], x[-1]], [D,D], label = 'D')
    plot(x, d, label = 'd')
    plot(x, L, label = 'L')
    xlabel('$x$ (mm)')
    legend()
### ---------------------------------------------------------------------------


#%% P vs t (Chamber pressure vs time)
###### Now, the second table on Nakka's spreadsheet: chamber pressure vs time
### Chamber cross sectional flow area
A_duct = (pi/4)*(dc**2 -(D**2-d**2))
### Propellant 'erosive burning' factor
G = G_star- A_duct/At
G = where(G>=0, G, 0.0)
### ideal throat area (m²) (should be equal to the throat area)
A_star0 = At0/1e6
A_star = At/1e6

M = propellant.M   # effective molecular mass of products (g/mol)
R = R_prime/M  # specific gas constant J/gK
k = propellant.k   # ratio of specific heats cP/cV, mixture
To = propellant.To
To_act = hc*To  # actual chamber temperature (K)
#c_star = 885.0  # characteritic exhaust velocity (m/s)
### I think this bellow is the ideal c_star. The definition is c* = p1 At/m_dot
c_star_ideal = sqrt(1000*R*To_act/k*((k+1)/2)**((k+1)/(k-1)))
### sound velocity at the throat 
vs_star = sqrt(1000*R*To_act*k)


N_P = int(1.2*N_web) #3000 # number of points to integrate burning rate versus pressure
### Chamber (stagnant) pressure Po (MPa)
Po = zeros(N_P)
Po[0] = Patm # Use MPa for burn rate equation
### burn rate (mm/s)
r = zeros_like(Po)
#r[0] = (1+kv*G[0])*a*Po[0]**n
r[0] = propellant.calculate_burning_rate(Po[0], G[0], kv)
### mass of propellant in the grain (kg)
mg = zeros_like(Po)
mg[0] = mg0/1000
# time since start of combustion
t = zeros_like(Po)
t[0] = 0.0
# volume of the grain (mm³)
Vg = zeros_like(Po)
Vg[0] = cyl_vol(D0, d0, L0) ### mm³
# Free volume in the chamber (m³)
Vfree = zeros_like(Po)
Vfree[0] = Va/1e6 -Vg[0]/1e9 ### m³ (Va is the chamber total volume, in cm³)
# mass generation rate of combustion products (kg/s)
dmdt_prod = zeros_like(Po)
dmdt_prod[0] = rho_p*r[0]*Ab0/1e6 #0.0
# mass flow rate through the nozzle (kg/s)
dmdt_flow = zeros_like(Po)
dmdt_flow[0] = 0.0
# mass storage rate in the chamber (kg/s)
dmdt_sto = zeros_like(Po)
dmdt_sto[0] = dmdt_prod[0]-dmdt_flow[0]
# mass of combustion products stored in the chamber (kg)
m_sto = zeros_like(Po)
m_sto[0] = 0.0
# ejection velocity (m/s) (column AI, in the spreadsheet)
#ve = zeros_like(Po)
#ve[0] = 0.0 # m/s

### density of combustion products in chamber (kg/m³)
rho_prod = zeros_like(Po)
rho_prod[0] = 0.0

for i in range(1, N_web):
    if inhibit_ext:
        Vg[i] = cyl_vol(D0, d[i], L[i]) ### mm³
    else:
        Vg[i] = cyl_vol(D[i], d[i], L[i]) ### mm³
    Vfree[i] = Va/1e6 -Vg[i]/1e9 ### m³ (Va is the chamber total volume, calculated previously in cm³)
    mg[i] = rho_p*Vg[i]/1e6 # kg
    dt = dx/r[i-1] # s
    t[i] = t[i-1] + dt # s
    dmdt_prod[i] = -(mg[i]-mg[i-1])/dt # kg/s
    #dmdt_prod[i] = rho_p*r[i-1]*Ab[i]/1e6
    ## R was given in J/g.K, to convert to SI units R_SI = 1000 R [J/kg.K]
    #dmdt_flow[i] = 1e6*(Po[i-1]-Patm)*A_star/sqrt(1000*R*To_act/k)*(2/(k+1))**((k+1)/(2*(k-1))) #rho_prod[i]*(At/1e6)*ve[i] # kg/s
    dmdt_flow[i] = 1e6*(Po[i-1]-Patm)*(At[i]/1e6)/sqrt(1000*R*To_act/k)*(2/(k+1))**((k+1)/(2*(k-1))) #rho_prod[i]*(At/1e6)*ve[i] # kg/s
    dmdt_sto[i] = dmdt_prod[i]-dmdt_flow[i]
    m_sto[i] = m_sto[i-1]+dt*dmdt_sto[i]
    rho_prod[i] = m_sto[i]/Vfree[i]
    Po[i] = rho_prod[i]*R*To_act/1000 +Patm
    ### burn rate
    #r[i] = (1+kv*G[i])*a*Po[i]**n ### use Po in MPa
    r[i] = propellant.calculate_burning_rate(Po[i], G[i], kv)

t_burn = t[i]
Pmax = max(Po)
print(f't_burn = {t_burn:.5g} s, Pmax = {Pmax:.3g} MPa')

for i in range(N_web, N_P):
    Vg[i] = 0.0
    Vfree[i] = Va/1e6
    mg[i] = 0.0
    #dt = dx/r[i-1] # s
    t[i] = t[i-1] + dt # s
    dmdt_prod[i] = 0.0
    #dmdt_flow[i] = 1e6*(Po[i-1]-Patm)*A_star/sqrt(1000*R*To_act/k)*(2/(k+1))**((k+1)/(2*(k-1))) #rho_prod[i]*(At/1e6)*ve[i] # kg/s
    dmdt_flow[i] = 1e6*(Po[i-1]-Patm)*(At[-1]/1e6)/sqrt(1000*R*To_act/k)*(2/(k+1))**((k+1)/(2*(k-1))) #rho_prod[i]*(At/1e6)*ve[i] # kg/s
    dmdt_sto[i] = -dmdt_flow[i]
    m_sto[i] = m_sto[i-1]+dt*dmdt_sto[i]
    rho_prod[i] = m_sto[i]/Vfree[i]
    Po[i] = rho_prod[i]*R*To_act/1000 +Patm
    ### burn rate
    r[i] = 0.0


#%%% plot P vs t and accessory variables
#### combustion variables: volumes, mass rates, burning rate, etc.
if plot_combustion_vs_t:
    figure('Integration variables', figsize = (10,8))
    # plot(t, Vg/1e9, label = '$V_g$ (m³)')
    # plot(t, Vfree, label = '$V_{free}$ (m³)')
    # plot(t, mg, label = '$m_g$ (kg)')
    # plot(t, dmdt_prod, label = '$\\dot m_{prod}$ (kg/s)')
    # plot(t, dmdt_flow, label = '$\\dot m_{flow}$ (kg/s)')
    # plot(t, m_sto, label = '$m_{stored}$ (kg)')
    # #plot(t, rho_prod, label = '$\\rho_{prod}$ (kg/m³)')
    # plot(t, r, label = '$r$ (mm/s)')
    # xlabel('time (s)')
    #legend()
    ax1 = subplot(4,1,1)
    #plot(t, Vg/1e9, label = '$V_g$ (m³)')
    #plot(t, Vfree, label = '$V_{free}$ (m³)')
    #ylabel('volumes (m³)')
    plot(t, Vg/1e6, label = '$V_g$ (L)')
    plot(t, Vfree*1e3, label = '$V_{free}$ (L)')
    ylabel('volumes (L)')
    legend()
    title(prop_name)
    tick_params('x', labelbottom=False)
    ax2 = subplot(4,1,2, sharex=ax1)
    plot(t, mg, label = '$m_g$ (kg)')
    plot(t, m_sto*100, label = '$m_{stored}$ (dag)')
    ylabel('masses (kg, dag)')
    legend()
    tick_params('x', labelbottom=False)
    ax3 = subplot(4,1,3, sharex=ax1)
    plot(t, dmdt_prod, label = '$\\dot m_{prod}$ (kg/s)')
    plot(t, dmdt_flow, label = '$\\dot m_{flow}$ (kg/s)')
    ylabel('mass rates (kg/s)')
    legend()
    tick_params('x', labelbottom=False)
    ax4 = subplot(4,1,4, sharex=ax1)
    #plot(t, rho_prod, label = '$\\rho_{prod}$ (kg/m³)')
    #plot(t, r, label = '$r$ (mm/s)')
    plot(t, rho_prod, label = '$\\rho_{prod}$ (kg/m³)')
    plot(t, r, label = '$r$ (mm/s)')
    ylabel("burn rate (mm/s)\n gas density (kg/m³)")
    legend()
    xlabel('time (s)')
    legend()
####

### Block to plot P vs time
#print(f"Using propellant {fuel_name}")
if plot_P_vs_t:
    figure('P_vs_t')
    plot(t, Po, label = prop_name)
    ylabel("$P_c$ (MPA)")
    xlabel("time (s)")
    legend()
###


#%%% First estimate of thrust vs time
### We can already make a first estimate of the thrust using
### F = ve*dmdt_flow
#F1 = c*dmdt_flow  
#### First approximation: 
# ve=zeros_like(t)
# ve[1:] = sqrt(2e6*(Po[1:]-Patm)/rho_prod[1:])
# ve[0] = ve[1]
# F1 = ve*dmdt_flow
# ### another estimate: F = Po*At
# #F2 = 1e6*(Po-Patm)*A_star
F2 = zeros_like(Po)
F2[:N_web] = (Po[:N_web]-Patm)*At
F2[N_web:] = (Po[N_web:]-Patm)*At[-1]
### Correct:
### F = m_dot v_2 +(P2-Patm)*A2
### v2 = sqrt( (2k/(k-1))R T1 [1-(P2/P1)^((k-1)/2)] +v1^2 )
### Let's assume P1 = Po, v1 = 0, T1 = To_act
ve = sqrt((2*k/(k-1))*1000*R*To_act*(1-(Patm/Po)**((k-1)/k)))
F1 = dmdt_flow*ve 

dmdt_flow_ave = mean(dmdt_flow)
P1 = mean(Po)
c_star = P1*At0/dmdt_flow_ave
print(f'c_star = {c_star:.4g} m/s, c_star ideal = {c_star_ideal:.4g} m/s.')

if plot_F_vs_t:
    #figure('F_vs_t_1')
    fig_force, ax_force = subplots(num = 'F_vs_t_1')
    plot(t, F1, label = '$F_1$')
    ylabel('$F$ (N)')
    xlabel('time (s)')
    plot(t, F2, label = '$F_2$')
    legend()

### calculating total impulse using trapezoidal rule (integration)
I = zeros_like(t)
F_max = max(F1)
### mean(F) does not correspond exactly to the average value, because t is
### not uniformly spaced. Calculating with and integral reveals that F_avg = It/t_thrust
pulse_started = False
for i in range(1,N_P):
    I[i] = (F1[i]+F1[i-1])*(t[i]-t[i-1])/2
    if F1[i] > 0.01*F_max:
        if not(pulse_started):
            t_start = t[i]
            pulse_started = True
        else:
            t_pulse = t[i]
It = sum(I)
t_pulse -= t_start
F_avg = It/t_pulse
print(f'Fmax = {F_max:.5g} N, Favg = {F_avg:.5g} N, t_thrust = {t_pulse:.5g} s')
print(f'It = {It:.5g} N.s, Isp = {1000*It/(g*mg0):.3g} s')
def motor_class(It):
    """ Classification ('letter') of the motor. Valid up to 'Z' (41,943 kNs to 83,886 kNs).

    pameters:
        It is the total impulse, in N.s (a numeric value, int or float).
    return:
        Class: a string indicating the motor class, from 'Micro' to 'Z'.

    Ref.: https://en.wikipedia.org/wiki/Model_rocket_motor_classification#Motor_impulse_by_class

    Author: Hugo L. D. de Souza Cavalcante."""

    if It < 0.15625:
        return "Micro"
    if It < 0.3125:
        return "1/8A"
    elif It<0.625:
        return "1/4A"
    elif It<1.26:
        return "1/2A"
    elif It>83886000:
        return "Larger than Z"
    else:
        return chr(int(ceil(log2(2*(It)/5)))+65)

classification = motor_class(It)
print(f'Motor classification: {classification} {round(F_avg)}-0')

# %%  Thrust vs time
# %%% Nozzle geometry
######### Nozzle geometry
### Nozzle convergence half-angle beta (degrees)
beta = 89 #70 #89.0 #35.0
beta_rad = beta*pi/180.0
### Nozzle divergence half-angle alpha (degrees)
alpha = 85 #12 #89.0 #12.0
alpha_rad = alpha*pi/180.0
### Nozzle exit diameter (mm)
De = dc # mm. If truncated, you can enter a value. Affects the parameter L_nozzle_d (divergence lenght)
De = Dt0+1.0
L_nozzle_c = ((dc-Dt0)/2)/tan(beta_rad)
L_nozzle_d = ((De-Dt0)/2)/tan(alpha_rad)
L_nozzle_t = 2.0 #25.0 # mm 
L_nozzle = L_nozzle_c + L_nozzle_d + L_nozzle_t ## mm 
#########

### block to draw nozzle geometry
if draw_nozzle:
    # figure('Nozzle geometry', figsize = (3,3.5))
    fig, ax = subplots()
    
    plot([-L_nozzle_t/2-(Dtf-Dt0)/tan(alpha_rad)/2, L_nozzle_t/2+(Dtf-Dt0)/tan(beta_rad)/2],
         [Dtf/2, Dtf/2], lw = 2, color = 'red')
    plot([-L_nozzle_t/2-(Dtf-Dt0)/tan(alpha_rad)/2, L_nozzle_t/2+(Dtf-Dt0)/tan(beta_rad)/2],
         [-Dtf/2, -Dtf/2], lw = 2, color = 'red')
    plot([-L_nozzle_d-L_nozzle_t/2, -L_nozzle_t/2, L_nozzle_t/2, L_nozzle_t/2+L_nozzle_c], 
         [De/2, Dt0/2, Dt0/2, dc/2], lw = 2, color = 'black')
    plot([-L_nozzle_d-L_nozzle_t/2, -L_nozzle_t/2, L_nozzle_t/2, L_nozzle_t/2+L_nozzle_c], 
         [-De/2, -Dt0/2, -Dt0/2, -dc/2], lw = 2, color = 'black')
    
    ax_x_lims = ax.get_xlim()
    ax_y_lims = ax.get_ylim()
    ax.set_box_aspect((ax_y_lims[1]-ax_y_lims[0])/(ax_x_lims[1]-ax_x_lims[0]))
    grid()
    arrow(0,0, -L_nozzle/5, 0, width = 0.8)
    #fig.tight_layout()


# %%% Third spreadsheet in Nakka's SRM_2023: "Performance"
### more constants
#etanoz = 0.85 # nozzle efficiency
#Ae0 = etanoz*At0 # nozzle exit area (initial) (mm²)
### local sonic velocity
#c = sqrt(k*R*T) also known as vs
#Me = # Mach number at nozzle exit

### State equation for ideal gas: PV = n_moles R_prime T
### or P = rho_prod R T   # R is the specific gas constant R = R_prime/M_molar
### Speed of sound in an ideal gas: c = sqrt(k P/rho_prod)
### k = ratio of specific heats cP/cV
### c = sqrt(k R T) 
### Mach number: M = v/c
### Differential equation for pressure along nozzle
### dP (1-M²) = rho_prod v² (dA/A)  # A is the (transversal) area of the duct (nozzle)
### Properties equations for isentropic flow:
### To/T = (1 + (k-1)/2 M²)
### Po/P = (1 + (k-1)/2 M²)^(k/(k-1))
### rho_o/rho = (1+(k-1)/2 M²)^(1/(k-1))
### dmdt_flow = (Po/sqrt(R To/k)) A M (1 + (k-1)/2 M²)^((k+1)/(2(1-k)))
### or dmdt_flow = (Po/sqrt(R To/k)) A M (1/(1 + (k-1)/2 M²))^((k+1)/(2(k-1)))
### A/A* = (1/M)( ( (k+1)/2 )/( 1+(k-1)/2 M² ) )^((k+1)/(2(1-k)))
### or A/A* = (1/M)( ( 1+(k-1)/2 M² )/( (k+1)/2 ) )^((k+1)/(2(k-1)))
### or A/A* = (1/M)( ( 1+(k-1)/2 M² )/( 1+(k-1)/2 ) )^((k+1)/(2(k-1)))
### Velocity at the nozzle exit
### ve = sqrt(2 To R (k/(k-1)) (1-(Pe/Po)^((k-1)/2)) ) # Pe is the exit pressure
### Pe = Patm, if
### A*/Ae = ((k+1)/2)^(1/(k-1)) (Pe/Po)^(1/k) sqrt(((k+1)/(k-1)) (1-(Pe/Po)^((k-1)/k)))
###
###



### Nozzle area expansion ratio (Ae/At)
exprat0 = (De/Dt0)**2
exprat = (De/Dt)**2
### Nozzle exit area Ae
#Ae = exprat*At
Ae = mean(exprat*At)
### Optimum nozzle area (mm²)  (using the equation for Pe = Patm)
At_opt = Ae*((k+1)/2)**(1/(k-1)) * (Patm/Po)**(1/k) * sqrt(((k+1)/(k-1))*(1-(Patm/Po)**((k-1)/k)))

At_opt_ave = mean(At_opt)
At_opt_max = max(At_opt)
At_opt_smart = mean(At_opt[N_web//15:N_web])
### Optimum nozzle diameter (mm)
De_opt = sqrt(Ae/At_opt_smart)*Dt0

### Exit velocity ve (m/s) (assuming Pe = Patm) (already calculated)
#ve = sqrt(2*To_act*1000*R*(k/(k-1))*(1-(Patm/Po)**((k-1)/2)))
### temperature at exit 
### cP = cV + R, k = cP/cV, so cP = R(k/(k-1))
### To = T + v²/(2cP)
Te = To_act - (ve**2)/2/R/k*(k-1)/1000
### local sound velocity at exit
vse = sqrt(k*1000*R*Te)
### Mach number at exit
Me = ve/vse #3.0 ### ?
### Exit Pressure (MPa)
Pe = Po/(1+(k-1)/2*Me**2)**(k/(k-1))
Pe = where(Pe>=Patm, Pe, Patm)

### Correcting the exit velocity with the true exit pressure 
#ve = sqrt(1000*R*To_act/k)*((k+1)/2)**((3-k)/(2*k-2))
v2 = sqrt(2*To_act*1000*R*(k/(k-1))*(1-(Pe/Po)**((k-1)/2)))

T2 = To_act - (v2**2)/2/R/k*(k-1)/1000
vs2 = sqrt(k*1000*R*T2)
M2 = v2/vs2 
P2 = Po/(1+(k-1)/2*M2**2)**(k/(k-1))
P2 = where(P2>=Patm, P2, Patm)
v3 = sqrt(2*To_act*1000*R*(k/(k-1))*(1-(P2/Po)**((k-1)/2)))

T3 = To_act - (v3**2)/2/R/k*(k-1)/1000
vs3 = sqrt(k*1000*R*T3)
M3 = v3/vs3 
P3 = Po/(1+(k-1)/2*M3**2)**(k/(k-1))
P3 = where(P3>=Patm, P3, Patm)

v4 = sqrt(2*To_act*1000*R*(k/(k-1))*(1-(P3/Po)**((k-1)/2)))


#c = where(dmdt_flow > 0.01, v4 - (P3-Patm)*Ae/dmdt_flow , 0.0) 
#c = where(c>=0.0, c, Patm)
c = ve - (Pe-Patm)*Ae/dmdt_flow_ave


### Recalculating the force
### F = m_dot v_2 +(P2-Patm)*A2
F3 = dmdt_flow*v2 +(Pe-Patm)*Ae
F4 = dmdt_flow*v3 +(P2-Patm)*Ae
F5 = dmdt_flow*v4 +(P3-Patm)*Ae


#figure('Force corrected')
#plot(t,F3)
ax_force.plot(t, F3, label = '$F_3$')
ax_force.plot(t, F4, label = '$F_4$')
ax_force.plot(t, F5, label = '$F_5$')
#xlabel('time (s)')
#ylabel('$F$ (N)')
ax_force.legend()



### something are wrong here (area expansion ratio not affecting performance)
### We need to calculate the pressure evolution along the nozzle, with?
### dP (1-M²) = rho_prod v² (dA/A)  
# N_nozzle = 10000
# dx_nozzle = L_nozzle_d/(N_nozzle-1)
# x_nozzle = zeros(N_nozzle)
# D_nozzle = zeros_like(x_nozzle)
# D_nozzle[0] = Dt[N_web//15]/1000 # or Dt0
# Px = zeros_like(x_nozzle)
# Px[0] = Po[N_web//15]*1e6
# rho = zeros_like(x_nozzle)
# rho[0] = Px[0]/R/To_act/1000
# A = zeros_like(x_nozzle)
# A[0] = At[N_web//15]/1e6
# v = zeros_like(x_nozzle)
# #v[0] = sqrt(2e6*(Po[1:]-Patm)/rho_prod[1:])
# v[0] = sqrt(1000*R*To_act*k) # M = 1 ?
# vs = zeros_like(x_nozzle)
# vs[0] = sqrt(1000*R*To_act*k) 
# T = zeros_like(x_nozzle)
# T[0] = To_act
# rho_v_A_const = rho[0]*v[0]*A[0]/1000
# Mach = 1.1
# for i in range(1,N_nozzle):
#     D_nozzle[i] = D_nozzle[i-1] +2*dx_nozzle*tan(alpha_rad)/1000
#     A[i] = (pi/4)*D_nozzle[i]**2
#     dA = A[i]-A[i-1]
#     T[i] = To_act/(1+((k-1)/2)*Mach**2)
#     vs[i] = sqrt(1000*R*T[i]*k) 
#     rho[i] = Px[i-1]/R/T[i]/1000
#     v[i] = rho_v_const*A[i-1]/(rho[i]*A[i])
#     Mach = v[i]/vs[i]
#     Px[i] = Px[i-1] + rho_v_A_const*v[i]*dA/A[i]/(1-Mach**2)/1000
#     #rho[i] = Px[i]/R/T[i-1]
#     x_nozzle[i] += dx_nozzle

#%%% plot variables associated with the nozzle exit
if plot_nozzle_vars:
    figure('Nozzle-exit variables', figsize = (10,8))
    ax1 = subplot(4,1,1)
    plot(t, ve, label = '$v_e$ (m/s)')
    plot(t, vse, label = '$v_s$ (m/s)')
    ylabel('exit and sound velocity (m/s)')
    legend()
    title(prop_name)
    tick_params('x', labelbottom=False)
    ax2 = subplot(4,1,2, sharex=ax1)
    plot(t, Te, label = '$T_e$ (K)')
    plot([t[0], t[-1]], [To, To], label = '$T_o$ (theoric) (K)')
    plot([t[0], t[-1]], [To_act, To_act], label = '$T_o$ actual (K)')
    ylabel('exit and chamber temperatures (K)')
    legend()
    tick_params('x', labelbottom=False)
    ax3 = subplot(4,1,3, sharex=ax1)
    plot(t, Me, label = 'Mach number')
    
    ylabel('Mach number')
    #ylabel('exit pressure (MPa) and Mach number')
    #legend()
    tick_params('x', labelbottom=False)
    ax4 = subplot(4,1,4, sharex=ax1)
    plot(t, Pe, label = '$P_e$ (MPa)')
    plot(t, P2, label = '$P_2$ (MPa)')
    plot(t, P3, label = '$P_3$ (MPa)')
    ylabel('exit pressure (MPa)')
    #legend()
    xlabel('time (s)')
    legend()
####



#%% Fim 
### It used to be necessary to plot figures in a noninteractive session
show()

