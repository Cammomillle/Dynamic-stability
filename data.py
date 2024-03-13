import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches

#************
#*** DATA ***
#************

# SI UNITS !!!
g = 9.81    # m/s^2
rho = 0.909 # at 10,000 feet 
V_0 = 28.809 # m/s

#******* Weights [kilograms * g] **********
W_w = 160*g     # wings
W_t_h = 10*g    # tail horizontal
W_t_v = 20*g    # tail vertical
W_fus = 90*g    # fuselage
W_motor = 20*g      # motor
W_batteries = 80*g  # batteries
W_gear = 10*g   # gears
W_flight_control = 30*g # commandes

#********* Fuselage ***********************
l_fus = 8.5     # length of the fuselage
width_fus = 0.6986   # Width of the fuselage
S_fuselage = 11.52   # Surface of the fuselage             A MODIFIER !!!!
z_nose = 0.334       # Height of the nose

#********** Wing **************************
S_w_total = 18      # total surface of the wing
S_w = S_w_total/2   # surface of half a wing    
b_w = 20            # Span of the wing
AR_w = b_w**2/S_w_total # AR of the wing
lambda_w = 0.4      # Taper ratio of the wing
c_w_root = 1.28571  # Chord of the wing at the root
c_w_tip = 0.51429   # Chord of the wing at the tip
c_mac_w = 2/3 * c_w_root * (1 + lambda_w + lambda_w**2)/(1 + lambda_w)  # mean aerodynamic chord of the wing
sweep_w = 0        # Sweep angle of the wing (LE sweep angle)
                                               
a0_w = 6.31206     # dcl_w/dalpha 
a_w = 5.62704      # dCL_w/dalpha 
c_m_w = -0.18      # Moment coefficient of the wing
z_wing = 0.85      # Height of the wing
aoa_w = 0          # AoA in configuration 1
L_w = 6864.655     # Lift from the wing in configuration 1
D_w = 198.99029    # Drag from the wing in configuration 1

def compute_lift_curve_slope(AR, cl_alpha, sweepback_angle_half, beta=1): # a for swept wwings/empennage !! other formula for unswept
    k=beta*cl_alpha/(2*np.pi)
    Lambda_beta=sweepback_angle_half
    return 2*np.pi/(beta*(2/(beta*AR)+np.sqrt(1/(k**2*(np.cos(Lambda_beta))**2)+(2/(beta*AR))**2)))

#*********** Fin (=empennage vertical) **************************
c_v_root = 1.17    # chord at the root of the fin
c_v_tip = 0.77     # chord at the tip of the fin
b_v = 1.5          # span of the fin
S_fin = (c_v_root+c_v_tip)*b_v/2   # surface of the fin (vertical tail)
lambda_v = c_v_tip/c_v_root         # taper ratio of the fin
c_mac_v = 2/3 * c_v_root * (1 + lambda_v + lambda_v**2)/(1 + lambda_v)  # mean aerodynamic chord of the fin
AR_v = b_v**2/S_fin   # AR of the horizontal empennage
sweep_v = np.arctan((0.75*c_v_root + np.tan(8*np.pi/180)*b_v - 0.75*c_v_tip)/b_v)  # Sweep angle of the fin (LE sweep angle)
sweep_v_half = np.arctan((0.5*c_v_root + np.tan(8*np.pi/180)*b_v - 0.5*c_v_tip)/b_v) # Sweep angle of the fin at the half chord

dclF_dbeta = 6.36    # slope of the CL curve for the fin in yaw motion (2D) = a0_v
dCLF_dbeta = compute_lift_curve_slope(AR_v*2, dclF_dbeta, sweep_v_half, beta=1)  # slope of the CL curve for the fin in yaw motion (3D) = a_v

#*********** Tail (=empennage horizontal) **************************
c_h_root = 0.55   # chord at the root of the horizontal empennage
c_h_tip = 0.35    # chord at the tip of the horizontal empennage
b_h = 4           # span of the horizontal empennage

S_h = (c_h_root+c_h_tip)*b_h/2  # surface of the horizontal empennage
lambda_h = c_h_tip/c_h_root     # taper ratio of the horizontal empennage
c_mac_h = 2/3 * c_h_root * (1 + lambda_h + lambda_h**2)/(1 + lambda_h)  # mean aerodynamic chord of the tail
AR_h = b_h**2/S_h   # AR of the horizontal empennage
sweep_h = np.arctan((0.75*(c_h_root-c_h_tip))/(b_h/2)) # Sweep angle of the horizontal empennage
sweep_h_half = np.arctan(0.5*(c_h_root-c_h_tip)/(b_h/2)) # Sweep angle of the horizontal empennage at the half chord

#Re=879.000 for the horizontal stabilizer. With airfool tools, at Re=1.000.000:
CD_0_h = 0.014          # same as the wing, from conceptual design, sld 61
dclh_alpha_h = 6.36     # slope of the CL curve for the tail (2D) for now, same as the fin !!!Re
dCLh_alpha_h = compute_lift_curve_slope(AR_h, dclh_alpha_h, sweep_h_half, beta=1)  # slope of the C_L of the tail with alpha of the tail   

z_tail = 1.0216 + b_v  # Height of the tail

#********** Position of CGs [meters] *******
# Compute the x-location of the aerodynamic center from the nose of the sailplane (formula for linearly-tapered wings)
def compute_x_ac(x_LE, b, taper, sweep, c_mac): 
    return x_LE + b/6*((1+2*taper)/(1+taper))*np.tan(sweep) + 0.25*c_mac

x_debut_wing = 2.95
x_w = compute_x_ac(x_debut_wing, b_w, lambda_w, sweep_w, c_mac_w) # Assumption : CG of the wing at the position of the AC of the wing

x_debut_tail = l_fus-c_v_root
x_t_v = compute_x_ac(x_debut_tail, b_v, lambda_v, sweep_v, c_mac_v) # Assumption : CG of the fin at the position of the AC of the fin

x_debut_tail_h = l_fus + np.tan(8*np.pi/180)*b_v - c_v_tip
x_t_h = compute_x_ac(x_debut_tail_h, b_h, lambda_h, sweep_h, c_mac_h) # Assumption : CG of the tail at the position of the AC of the tail

x_crew1 = 1.435       # crew member 1
x_crew2 = 2.575       # crew member 2

x_b_w = x_crew1     # ballast au CG du crew 1
x_b_t = x_t_v       # ballast au CG du fin

x_fus = 0.45*l_fus         # fuselage
x_motor = 4.               # motor                              A MODIFIER !!!!
x_batteries = 4.4          # batteries                        A MODIFIER !!!!
x_main_gear = 3.5          # main gear                          A MODIFIER !!!!
x_gear_2 = l_fus-0.55      # last gear                          A MODIFIER !!!!
x_gear = (x_main_gear*0.75+x_gear_2*0.25)/(x_main_gear+x_gear_2)
x_flight_control = x_crew1-0.2 # assumption of 20cm between first pilots and commands

#********** Pitch stability **************
K_n_low = 0.05      # Stability min to be stable
K_n_high = 0.25     # Stability max to not be too stable

#********** Yaw stability ****************
hf1 = 0.702+0.285      # see Conceptual design slide 57
hf2 = 0.2356           
bf1 = 0.3346+0.342
bf2 = 0.2
h_f_max = 0.974
