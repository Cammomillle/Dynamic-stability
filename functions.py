from data import *
from scipy.optimize import curve_fit, fmin
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
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w_le) # x-location of the LE of the projected c_mac of the wing
    h_0 = (x_ac_w-x_le_wing)/c_mac_w
    a1 = dCLh_alpha_h
    l_T = x_ac_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)
    m = 2*(z_tail-z_wing)/b_w
    grad_downwash = 1.75*a_w/(np.pi * AR_w * ((2*lambda_w*l_T)/b_w)**0.25*(1+abs(m)))
    m_fus = x_ac_w/l_fus
    k_fus = k_fus_table(m_fus)
    grad_Cmfus_Clw = k_fus*width_fus**2*l_fus/(S_w*c_mac_w*a_w)
    h_n = h_0+V_t_bar*a1/a_w*(1-grad_downwash) - grad_Cmfus_Clw
    return h_n

def compute_h(x_cg):
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w_le) # x-location of the LE of the projected airfoil of the wing
    h=(x_cg-x_le_wing)/c_mac_w # by definition !
    return h

def CL_alpha_wing(M):
    beta=np.sqrt(1-M**2)
    k=a0_w/(2*np.pi)
    CL_alphaW=(2*np.pi*AR_w)/(2+np.sqrt((AR_w*beta/k)**2*(1+np.tan(sweep_w_half)**2/beta**2)+4))
    return CL_alphaW

def downwash():
    K_A = 1/AR_w - 1/(1+AR_w**1.7)
    K_lambda = (10-3*lambda_w)/7
    K_H = (1-h_H/b_w)/(np.sqrt(2*l_H/b_w))
    grad_down_M0 = 4.44*(K_A*K_lambda*K_H*np.cos(sweep_w))**1.19
    grad_down = grad_down_M0*1      # approx that CL_alpha,wM/CL_alpha,wM=0 ~ 1
    return grad_down 

def body_sum(M, x_cg):  # computation of the data required to compute the body effect on downwash for CM_alpha
    dXi = 0.1
    frames = np.arange(pe.pos[0]/1000, pe.pos[-1]/1000, dXi)

    def hyperbole(x, a, b, c):
        return b/x**a + c

    figure11 = np.genfromtxt('LongData/Figure1_1.csv', delimiter=',')
    Rfigure11 = curve_fit(hyperbole, figure11[:,0], figure11[:,1])[0]
    figure12 = np.genfromtxt('LongData/Figure1_2.csv', delimiter=',')
    Rfigure12 = curve_fit(hyperbole, figure12[:,0], figure12[:,1])[0]
    
    global_dwash = downwash()

    ret = 0.0
    for frame in frames:
        Wf = 2*pe.b_i(frame*1000)/1000

        if frame < x_debut_wing:
            dwash = (1 - global_dwash)*(x_debut_wing - frame)/l_H
        elif frame > x_debut_wing + c_w_root:
            xicf = (frame - x_debut_wing - c_w_root)/c_w_root
            dwash = hyperbole(xicf, *Rfigure11) * CL_alpha_wing(0) / (0.08 * 180 / np.pi)
        else:
            dwash = 0.0

        ret += Wf**2 * dwash * dXi
    
    return ret

def compute_K_n_bis(W_b,x_b,W_crew1,W_crew2):
    d = 0.4956851278
    eta_h = 0.9
    K_WB = 1-0.25*(d/b_w)**2+0.025*(d/b_w)
    a_w = CL_alpha_wing(0)
    grad_downwash = downwash()
    x_cg = compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    CL_alpha_WB = K_WB*a_w 
    M_alpha = 0.5*rho*V0**2 / 36.5 * body_sum(0, x_cg)
    d_xac = -M_alpha / (0.5*rho*V0**2 * S_w_total * CL_alpha_wing(0))
    x_acwb = x_ac_w+d_xac
    num = x_acwb + dCLh_alpha_h/CL_alpha_WB * eta_h * S_h/S_w_total * x_ac_h * (1-grad_downwash)
    den = 1+dCLh_alpha_h/CL_alpha_WB * eta_h * S_h/S_w_total * (1-grad_downwash)
    x_ac = num/den
    Kn = -(x_cg-x_ac)/c_mac_w
    
    if(Kn >= K_n_low and Kn <= K_n_high):
        is_kn_ok = True
    else:
        is_kn_ok = False

    return Kn, is_kn_ok


