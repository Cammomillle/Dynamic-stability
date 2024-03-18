from data import *

#***************************************************************************
#*** Functions to compute the stability and equilibrium of the sailplane ***
#***************************************************************************

# Total mass of the sailplane
def compute_total_mass(W_b, W_crew1, W_crew2):
    return W_b + W_w + W_crew1 + W_crew2 + W_t_h + W_t_v + W_fus + W_motor + W_batteries + W_gear + W_flight_control

# Total mass of the empty sailplane (i.e. without crew and ballast)
def compute_total_mass_empty():
    return W_w + W_t_h + W_t_v + W_fus + W_motor + W_batteries + W_gear + W_flight_control

# Compute the center of gravity from the nose 
def compute_x_cg(W_b, x_b, W_crew1, W_crew2): 
    W_tot = compute_total_mass(W_b, W_crew1, W_crew2)
    return (W_b*x_b + W_crew1*x_crew1 + W_crew2*x_crew2 + W_w*x_w + W_t_h*x_t_h + W_t_v*x_t_v + W_fus*x_fus + W_motor*x_motor + W_batteries*x_batteries + W_gear*x_gear + W_flight_control*x_flight_control)/W_tot

# Compute the center of gravity from the nose for empty configuration (i.e. without crew and ballast)
def compute_x_cg_empty():
    W_tot = compute_total_mass_empty() 
    return (W_w*x_w + W_t_h*x_t_h + W_t_v*x_t_v + W_fus*x_fus + W_motor*x_motor + W_batteries*x_batteries + W_gear*x_gear + W_flight_control*x_flight_control)/W_tot

# Compute the x-location of the aerodynamic center from the nose of the sailplane (formula for linearly-tapered wings)
def compute_x_ac(x_debut, b, taper, sweep, c_mac): 
    return x_debut + b/6*((1+2*taper)/(1+taper))*np.tan(sweep) + 0.25*c_mac

# Compute the x-location of the projection of the c_mac on the fuselage from the nose of the sailplane (formula for linearly-tapered wings)
def compute_x_mac(x_debut, b, taper, sweep): 
    return x_debut + b/6*((1+2*taper)/(1+taper))*np.tan(sweep)
    
def k_fus_table(m_fus): # interpolation of table from slide 46 in Conceptual design course
    return np.exp(m_fus)

# Stability margin, for t-tail must be ]0.1;0.2[, can be a bit relaxed
def compute_K_n(h_n, h):
    K_n = h_n - h
    if(K_n >= K_n_low and K_n <= K_n_high):
        is_kn_ok = True
    else:
        is_kn_ok = False
    return K_n, is_kn_ok

def compute_hn(x_cg):  # hn is computed to define pitch stability, see Conceptual design slides 47-49
    x_ac_w = compute_x_ac(x_debut_wing, b_w, lambda_w, sweep_w, c_mac_w)
    x_ac_t_h = compute_x_ac(x_debut_tail_h, b_h, lambda_h, sweep_h, c_mac_h) 
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w) # x-location of the LE of the projected c_mac of the wing
    x_le_tail = compute_x_mac(x_debut_tail_h, b_h, lambda_h, sweep_h) # x-location of the LE of the projected c_mac of the tail
    
    h_0 = (x_ac_w-x_le_wing)/c_mac_w
    a1 = dCLh_alpha_h
    l_T = x_ac_t_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)
    m = 2*(z_tail-z_wing)/b_w
    grad_downwash = 1.75*a_w/(np.pi * AR_w * ((2*lambda_w*l_T)/b_w)**0.25*(1+abs(m)))
    m_fus = x_ac_w/l_fus
    k_fus = k_fus_table(m_fus)
    grad_Cmfus_Clw = k_fus*width_fus**2*l_fus/(S_w*c_mac_w*a_w)
    h_n = h_0+V_t_bar*a1/a_w*(1-grad_downwash) - grad_Cmfus_Clw
    return h_n

def compute_h(x_cg):
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w) # x-location of the LE of the projected airfoil of the wing
    h=(x_cg-x_le_wing)/c_mac_w # by definition !
    return h


print(compute_total_mass(15*g, 123*g, 123*g)/9.81)