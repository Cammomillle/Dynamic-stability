import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys
import pos_ellipse as pe
import matplotlib.patches as mpatches

#************
#*** DATA ***
#************

# SI UNITS !!!
g = 9.81    # m/s^2
rho = 0.909 # at 10,000 feet 
V0 = 28.809 # m/s
gamma = 1.4 

#******* Weights [kilograms * g] **********
W_w = 160*g     # wings
W_t_h = 10*g    # tail horizontale
W_t_v = 20*g    # tail verticale
W_fus = 90*g    # fuselage
W_motor = 20*g      # motor
W_batteries = 80*g  # batteries
W_gear = 10*g   # gears
W_flight_control = 30*g # commandes

#********* Fuselage ***********************
l_fus = 8.5     # length of the fuselage
width_fus = 0.6986   # Width of the fuselage
S_fuselage = 12.0816 # Surface of the fuselage             A MODIFIER !!!!
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
sweep_w = -1.105*np.pi/80       # Sweep angle of the wing (quarter chord)
sweep_w_half = -2.2089*np.pi/180 # Sweep angle of the wing (half chord)
sweep_w_le = 0  # sweep angle of the wing at the le
dihedral_w = 1*np.pi/180 
theta_w = 1*np.pi/180

a0_w = 6.31206     # dcl_w/dalpha 
a_w = 5.62704      # dCL_w/dalpha 
alpha_0 = -5.88608*np.pi/180 # zero lift angle [radians]
c_m_w = -0.18      # Moment coefficient of the wing
z_wing = 0.85      # Height of the wing
aoa_w = 0          # AoA in configuration 1
e_w = 0.9335        # Ostwald number
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
y_mac_v = 0.698    # à vérifier !!!
AR_v = b_v**2/S_fin   # AR of the horizontal empennage
sweep_v = np.arctan((0.75*c_v_root + np.tan(8*np.pi/180)*b_v - 0.75*c_v_tip)/b_v)  
sweep_v_half = np.arctan((0.5*c_v_root + np.tan(8*np.pi/180)*b_v - 0.5*c_v_tip)/b_v) # Sweep angle of the fin at the half chord
sweep_v_le = 22.16*np.pi/180 # Sweep angle of the fin (LE sweep angle)

dclF_dbeta = 6.36    # slope of the CL curve for the fin in yaw motion (2D) = a0_v
dCLF_dbeta = compute_lift_curve_slope(AR_v*2, dclF_dbeta, sweep_v_half, beta=1)  # slope of the CL curve for the fin in yaw motion (3D) = a_v

#**** Rudder (fin control surface) ****
b_r = 1.5
c_r_root = 0.33
c_r_tip = 0.21
c_ratio_root = c_r_root/c_v_root
c_ratio_tip = c_r_tip/c_v_tip
c_ratio = (c_ratio_root + c_ratio_tip)/2 # assumption: mean taken from graph
S_r = (c_r_root+c_r_tip)*b_r/2
AR_r = b_r**2/S_r
sweep_r_half = np.arctan((c_r_root/2+c_r_tip/2)/b_r)

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
sweep_h_le = 5.71*np.pi/180
theta_h = -1.45*np.pi/180

#Re=879.000 for the horizontal stabilizer. With airfool tools, at Re=1.000.000:
CD_0_h = 0.014          # same as the wing, from conceptual design, sld 61
dclh_alpha_h = 6.36     # slope of the CL curve for the tail (2D) for now, same as the fin !!!Re
dCLh_alpha_h = compute_lift_curve_slope(AR_h, dclh_alpha_h, sweep_h_half, beta=1)  # slope of the C_L of the tail with alpha of the tail   

i_h = -2*np.pi/180      # horizontal stabilizer incidence
z_tail = 1.0216 + b_v   # Height of the tail

#********** Position of CGs [meters] *******
# Compute the x-location of the aerodynamic center from the nose of the sailplane (formula for linearly-tapered wings)
def compute_x_ac(x_LE, b, taper, sweep, c_mac): 
    return x_LE + b/6*((1+2*taper)/(1+taper))*np.tan(sweep) + 0.25*c_mac

x_debut_wing = 2.9
x_ac_w = compute_x_ac(x_debut_wing, b_w, lambda_w, sweep_w_le, c_mac_w)
x_w = x_debut_wing + 0.428
theta_w = 1.4/180*np.pi
y_ac_root_w = -np.tan(theta_w)*c_w_root/4 + 0.65

x_debut_tail = l_fus-c_v_root
x_ac_v = compute_x_ac(x_debut_tail, b_v, lambda_v, sweep_v_le, c_mac_v)
x_t_v = x_debut_tail + 0.62

x_debut_tail_h = l_fus + np.tan(8*np.pi/180)*b_v - c_v_tip
x_ac_h = compute_x_ac(x_debut_tail_h, b_h, lambda_h, sweep_h_le, c_mac_h)
x_t_h = x_debut_tail_h + 0.276
x_crew1 = 1.435       # crew member 1
x_crew2 = 2.575       # crew member 2

x_b_w = x_crew1     # ballast au CG du crew 1
x_b_t = x_t_v       # ballast au CG du fin

x_fus = 3.3                # fuselage
x_motor = 4.               # motor  
x_batteries = 4.8          # batteries     
x_main_gear = 3.5          # main gear                          
x_gear_2 = l_fus-0.55      # last gear                          
x_gear = (x_main_gear*0.75+x_gear_2*0.25)/(x_main_gear+x_gear_2)
x_flight_control = x_crew1-0.2 # assumption of 20cm between first pilots and commands

#********** Aerodynamics **************
S_N = 0.5308287436110252 # maximum section area

#********** Pitch stability **************
K_n_low = 0.05      # Stability min to be stable
K_n_high = 0.25     # Stability max to not be too stable

#********** Yaw stability ****************
hf1 = 0.702+0.285      # see Conceptual design slide 57
hf2 = 0.2356           
bf1 = 0.3346+0.342
bf2 = 0.2
h_f_max = 0.974

# ---------- Dynamic stability ------------
alpha_e = 0     # equilibrium aoa
D_e = 170.378 #N

h_H = 1.646     # vertical distance between x_ac of tail and wing   !!! à vérifier !!!
l_H = x_ac_h - x_ac_w # horitontal distance between x_ac of tail and wing
z_V = 0.2      # A MODIFIER !!!!!
z_w = pe.ml_i(x_ac_w*1000)/1000 - y_ac_root_w # vertical distance positive downward between quarter chord of wing and sailplane centerline -> A MODIFIER ???
d = 0.4956851278    # = sqrt(averag fusegale cross sectional area/0.7854) from slide 13 lat. deriv. !! value that could change !!! 
d_fus = 0.637       # fuselage diameter (mean value)
alpha_dcl = -0.63   # fin flap value (extracted from slide 65)
alpha_dcl_ratio = 1.05 # fin flap value ratio (extracted from slide 65; conservative value)

eta = 20*np.pi/180  # max elevator deflection   
df = 30*np.pi/180   # max rudder deflection                               

Ix = 3383.31944 + 18.580411
Iy = 1320.1079 + 1331.37325
Iz = 4636.2118 + 1319.354609
Ixz = 137.65305


