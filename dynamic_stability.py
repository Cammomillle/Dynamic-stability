import numpy as np 
from matplotlib import pyplot as plt 
import control as ct 
import pos_ellipse as pe
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from data import *
from functions import *

"""
THE FOLLOWING VALUES MUST CHANGE !!
"""
CD0=0.014
T=-4.8+273.15 #temperature at 10 000 ft ! ->    shit ass guess !!!!
nu = 1.694*10**(-5) #Ns/m^2
#rho must be a function of altitude 
#z_F in data (dynamic stability part)
# a_rudder to verify
# "S" in derivatives commputation, not sure wether it is S_w_total or S_w

"""  END OF SUSPECT VALUES """

W_b_w = 0*g     # ballasts at wing
W_b_t = 0*g     # ballasts at tail
W_b = W_b_w + W_b_t  # ballasts total weight

W_crew1 = 80*g   # 1 st crew
W_crew2 = 80*g   # 2 nd crew

#******* CG of the ballasts **********
x_b = 0
if(W_b!=0):
    x_b = (W_b_w*x_b_w + W_b_t*x_b_t)/W_b

R=287 #J/kg/k pour l'air
W=compute_total_mass(W_b,W_crew1,W_crew2)
eta_h = 0.9 # slide 19 DATCOM (entre 0.9-1)

"-------------------------------------------------Utilities--------------------------------------"
def CL_alpha_wing(M):
    beta=np.sqrt(1-M**2)
    k=a0_w/(2*np.pi)
    CL_alphaW=(2*np.pi*AR_w)/(2+np.sqrt((AR_w*beta/k)**2*(1+np.tan(sweep_w_half)**2/beta**2)+4))
    return CL_alphaW

def CL_alpha_tail_horizontal(M):
    beta=np.sqrt(1-M**2)
    k=dclh_alpha_h/(2*np.pi)
    CL_alphaH=(2*np.pi*AR_h)/(2+np.sqrt((AR_h*beta/k)**2*(1+np.tan(sweep_h_half)**2/beta**2)+4))
    return CL_alphaH

def CL_beta_tail_vertical(M):
    beta=np.sqrt(1-M**2)
    k=dclF_dbeta/(2*np.pi)
    CL_alphaV=(2*np.pi*AR_v*2)/(2+np.sqrt((AR_v*2*beta/k)**2*(1+np.tan(sweep_v_half)**2/beta**2)+4))    # we multiply the AR_v by two because the formula is valide for a full wing
    return CL_alphaV/2

def CL_beta_rudder(M):
    beta=np.sqrt(1-M**2)
    k=dclF_dbeta/(2*np.pi)
    CL_alphaR=(2*np.pi*AR_r*2)/(2+np.sqrt((AR_r*2*beta/k)**2*(1+np.tan(sweep_r_half)**2/beta**2)+4))    # we multiply the AR_v by two because the formula is valide for a full wing
    return CL_alphaR/2

def compute_CL(V_0):
    return 2*W/(rho*V_0**2*S_w_total)

def compute_response(A,B,C,D,input_,T,X0):
    sys=ct.StateSpace(A,B,C,D)
    response=ct.forced_response(sys,T,input_,X0)
    
"------------------------------------------------Matrices----------------------------------------"
def compute_lat_matrices(m,Ix,Iy,Ixz,Iz,Y_p,Y_r,Y_v,Ue,We,theta_e,L_p,L_r,L_v,N_p,N_r,N_v,Y_ksi,Y_zeta,L_ksi,L_zeta,N_ksi,N_zeta):
    A1=np.matrix([[m,0,0,0,0],
                  [0,Ix,-Ixz,0,0],
                  [0,-Ixz,Iz,0,0],
                  [0,0,0,1,0],
                  [0,0,0,0,1]])
    
    B1=np.matrix([[-Y_v,-(Y_p+m*We),-(Y_r-m*Ue),-m*g*np.cos(theta_e),-m*g*np.sin(theta_e)],
                 [-L_v,-L_p,-L_r,0,0],
                 [-N_v,-N_p,-N_r,0,0],
                 [0,-1,0,0,0],
                 [0,0,-1,0,0]])
    
    C1=np.matrix([[Y_ksi,Y_zeta],
                 [L_ksi,L_zeta],
                 [N_ksi,N_zeta],
                 [0,0],
                 [0,0]])

    A_1=-np.linalg.inv(A1) @ B1
    B_1=np.linalg.inv(A1) @ C1
    return A_1, B_1

def compute_long_matrices(m,theta_e,Iy,X_w_dot,Z_w_dot,M_w_dot,X_u,X_w,X_q,Z_u,Z_w,Z_q,M_u,M_w,M_q,X_eta,X_tau,Z_eta,Z_tau,M_eta,M_tau,U_e,W_e):
    A2=np.matrix([[m, -X_w_dot, 0, 0],
                 [0, (m-Z_w_dot), 0, 0],
                 [0, -M_w_dot, Iy, 0],
                 [0, 0, 0, 1]])
    
    B2=np.matrix([[-X_u, -X_w, -(X_q-m*W_e), m*g*np.cos(theta_e)],
                 [-Z_u, -Z_w, -(Z_q+m*U_e), m*g*np.sin(theta_e)],
                 [-M_u, -M_w, -M_q, 0], 
                 [0, 0, -1, 0]])
    
    C2=np.matrix([[X_eta, X_tau], 
                  [Z_eta, Z_tau],
                  [M_eta, M_tau], 
                  [0, 0]])

    A_2=-np.linalg.inv(A2) @ B2
    B_2=np.linalg.inv(A2) @ C2

    return A_2, B_2

def long_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha_e):
    CD_u, CL_u, CM_u = u_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha_e,alpha_0)
    CD_alpha, CL_alpha, CM_alpha = alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0)
    CD_q, CL_q, CM_q = q_derivatives(W_b,x_b,W_crew1,W_crew2,V0)
    CD_alpha_dot, CL_alpha_dot, CM_alpha_dot = alpha_dot_derivatives(W_b,x_b,W_crew1,W_crew2,V0)

    L_e = compute_total_mass(W_b, W_crew1, W_crew2)
    D_e = 170.378 #N
    C_ze = (-L_e*np.cos(alpha_e) - D_e*np.sin(alpha_e))/(0.5*rho*V0**2*S_w_total)
    C_xe = (L_e*np.sin(alpha_e) - D_e*np.cos(alpha_e))/(0.5*rho*V0**2*S_w_total)
    
    X_u = CL_u*np.sin(alpha_e) - CD_u*np.cos(alpha_e)
    X_w = 1/(np.cos(alpha_e))*(-C_ze+CL_alpha*np.sin(alpha_e) - CD_alpha*np.cos(alpha_e))
    X_q = CL_q*np.sin(alpha_e) - CD_q*np.cos(alpha_e)
    X_w_dot = 1/np.cos(alpha_e)*(CL_alpha_dot*np.sin(alpha_e)-CD_alpha_dot)

    Z_u = -(CL_u*np.cos(alpha_e) + CD_u*np.sin(alpha_e))
    Z_w = 1/np.cos(alpha_e)*(C_xe - CL_alpha*np.cos(alpha_e) - CD_alpha*np.sin(alpha_e))
    Z_q = -(CL_q*np.cos(alpha_e) + CD_q*np.sin(alpha_e))
    Z_w_alpha_dot = -1/np.cos(alpha_e)*(CL_alpha_dot*np.cos(alpha_e) + CD_alpha_dot*np.sin(alpha_e))

    M_u = CM_u
    M_q = CM_q
    M_w = 1/np.cos(alpha_e)*CM_alpha
    M_w_dot = 1/np.cos(alpha_e)*CM_alpha_dot

    return X_u, X_w, X_q, X_w_dot, Z_u, Z_w, Z_q, Z_w_alpha_dot, M_u, M_q, M_w, M_w_dot

def lat_derivatives():
    return
"----------------------------------------------Inertias------------------------------------------"    
def compute_inertias():
    return 

"----------------------------------------------DATCOM longitudinal derivatives--------------------------------"
def alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    " SLIDE 19 DATCOM "
    " DIAMETRE FUSELAGE PUE DU CUL ON A PRIS EPAISSEUR FUSELAGE A LA PLACE MAIS CEST NUL!!!!"
    
    dCLh_alpha_h = CL_alpha_tail_horizontal(0)
    a_w = CL_alpha_wing(0)
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)

    #---- CL derivative ----
    CL=compute_CL(V0)
    d=0.4956851278
    K_WB=1-0.25*(d/b_w)**2+0.025*(d/b_w)
    grad_downwash = downwash()
    CL_alpha_WB=K_WB*a_w #KWB*Cl_alpha_wing
    CL_alpha=CL_alpha_WB+dCLh_alpha_h*eta_h*S_h/S_w_total*(1-grad_downwash)
    M_alpha = 0.5*rho*V0**2 / 36.5 * body_sum(0, x_cg)
    d_xac = -M_alpha / (0.5*rho*V0**2 * S_w_total * CL_alpha_wing(0))
    x_acwb=x_ac_w+d_xac
    num=x_acwb+dCLh_alpha_h/CL_alpha_WB * eta_h * S_h/S_w_total * x_ac_h * (1-grad_downwash)
    den=1+dCLh_alpha_h/CL_alpha_WB * eta_h * S_h/S_w_total * (1-grad_downwash)
    x_ac=num/den
    
    #--- CM derivative ----
    K = -(x_cg-x_ac)/c_mac_w
    #print('Stability margin', K)
    CM_alpha=-K*CL_alpha

    #--- CD derivative ----
    CD_alpha=2*CL*CL_alpha/(np.pi*AR_w*e_w)

    return CD_alpha, CL_alpha, CM_alpha


def u_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha,alpha_0):
    a=np.sqrt(R*gamma*T)
    Mach=V0/a 

    #---- CL derivative ----
    CL_alpha_for_mach,CD_alpha,CM_alpha=alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0) #CL_alpha for a given mach M !
       
    CL_for_mach=CL_alpha_for_mach*(alpha-alpha_0)
    CL_u=Mach**2/(1-Mach**2) * CL_for_mach # ATTENTION recompute CL for the given mach !!
    
    V0_plus=V0+2 #an increment on the speed
    Mach_plus=V0_plus/a
    dM=Mach_plus-Mach
    CL_alpha_plus,_,_=alpha_derivatives(W_b, x_b, W_crew1, W_crew2, V0_plus)
    CL_plus=CL_alpha_plus*(alpha-alpha_0)
    dCD_dM=(CL_plus**2-CL_for_mach**2)/(dM*np.pi*AR_w*e_w) #We assumed that CD0 was independent of M !! -> is it okay ?
    "dCD_dm doit être calculé comme suit : 1) calculer CL pour M et pour M+delta M et en tirer les CD associés depuis la drag polar puis déf de la dérivée "
    
    #---- CD derivative ----
    CD_u=Mach*dCD_dM
    "Here we will assume that x_ac_w doesn't move with the mach, which should be fine for low Mach n# (don't know any better)"

    #---- CM derivative ----
    x_ac_plus=x_w
    x_ac_for_mach=x_w 
    dxac_w_dM=x_ac_plus-x_ac_for_mach/dM
    CM_u=-CL_for_mach * dxac_w_dM

    return CD_u, CL_u, CM_u

def q_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    a=np.sqrt(R*gamma*T)
    M=V0/a #mach number
    
    #---- CD derivative ----
    CD_q=0 #DATCOM recommends to neglect it
    
    #---- CL derivative ----
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w) # x-location of the LE of the projected c_mac of the wing
    x_ac_w = x_w
    h_0 = (x_ac_w-x_le_wing)/c_mac_w
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    h=compute_h(x_cg)
    B=np.sqrt(1-(M*np.cos(sweep_w)**2))
    CL_qW=(AR_w+2*np.cos(sweep_w))/(AR_w*B+2*np.cos(sweep_w))*(1/2+2*(h-h_0))*CL_alpha_wing(M=0)
    l_T = x_ac_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)

    CL_qH=2*CL_alpha_tail_horizontal(M)*eta_h*V_t_bar
    CL_q=CL_qW+CL_qH
    
    #---- CM derivative ----
    num = (AR_w**3*(np.tan(sweep_w))**2)/(AR_w*B+6*np.cos(sweep_w)) + 3/B   # B est une correction de compressibilité => B ~ 1
    den = (AR_w**3*(np.tan(sweep_w))**2)/(AR_w+6*np.cos(sweep_w)) + 3
    num1 = AR_w*(2*((h-h_0)**2)+0.5*(h-h_0))/(AR_w+2*np.cos(sweep_w))
    num2 = 1/24*AR_w**3*(np.tan(sweep_w))**2/(AR_w+6*np.cos(sweep_w)) 
    num3 = 1/8
    K = 0.9 #voir graphe slide 42 DATCOM
    CM_qw0 = - K*a0_w*np.cos(sweep_w)*(num1 + num2 + num3)
    CM_qw = CM_qw0*(num/den)
    C_Lalphah = CL_alpha_tail_horizontal(M)
    CM_qh = -2*C_Lalphah*h*V_t_bar*l_T/c_mac_h  # slide 41 DATCOM c_mac -> c_mac_h ou c_mac_w ??
    CM_q = CM_qw + CM_qh

    return CD_q, CL_q, CM_q

def alpha_dot_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    a=np.sqrt(R*gamma*T)
    M=V0/a #mach number

    #---- CD derivative ----
    CD_alpha_dot = 0 # slide 43 DATCOM : negligible
    beta = np.sqrt(1-M**2)
    f_BA = 0.14/4 * beta*AR_w
    CL_g = f_BA*(np.pi/2*AR_w)/(-beta**2)

    #---- CL derivative ----
    CL_alpha_dot_w = 1.5*x_w/c_w_root*CL_alpha_wing(M) + 3*CL_g     # ici on peut négliger CL_alpha_dot_w : CL_alpha_dot = CL_alpha_dot_h
    x_cg = compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    l_T = x_ac_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)
    grad_downwash = downwash()
    CL_alpha_dot_h = 2*CL_alpha_tail_horizontal(M)*eta_h*V_t_bar*grad_downwash
    CL_alpha_dot = CL_alpha_dot_w + CL_alpha_dot_h

    #---- CM derivative ----
    CM_alpha_dot_w = 0
    CM_alpha_dot_h = -2*CL_alpha_tail_horizontal(M)*eta_h*V_t_bar*grad_downwash*l_T/c_mac_h # slide 47 DATCOM don't know if c_mac is c_mac_w or c_mac_h
    CM_alpha_dot = CM_alpha_dot_w + CM_alpha_dot_h

    return CD_alpha_dot, CL_alpha_dot, CM_alpha_dot

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

long_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha_e)

"----------------------------------------------DATCOM lateral derivatives--------------------------------"
def sidewash_beta():    # sidewash derivative wrt beta
    sigma_beta = -0.276 + 3.06*S_fin/S_w_total*1/(1+np.cos(sweep_v)) + 0.4*Z_w/d + 0.009*AR_v # slide 13 L13
    return sigma_beta

def Y_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    a_v = CL_beta_tail_vertical(M)
    a_r = CL_beta_rudder(M)

    # Y_v
    grad_sidewash = sidewash_beta()
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_F = x_ac_v - x_cg     
    Y_v_fin = -S_fin/S_w*a_v*(1-grad_sidewash)        # S = surface totale ou surface d'une demi aile ??!! 
    Y_v_wing = -0.0001*np.abs(dihedral_w)
    #Y_v_body = -2*Ki*                                # à faire !!!
    Y_v = Y_v_fin + Y_v_wing

    # Y_p            
    Y_p_fin = -2*S_fin/S_w*a_v*(1-grad_sidewash)*(z_F*np.cos(alpha_e)-l_F*np.sin(alpha_e))/b_v
    Y_p_wing = 0 # assumed negligible
    Y_p_body = 0
    Y_p = Y_p_fin + Y_p_wing + Y_p_body

    # Y_r
    Y_r_fin = 2*S_fin/S_w*a_v*(1-grad_sidewash)*(z_F*np.sin(alpha_e)+l_F*np.cos(alpha_e))/b_v
    Y_r_body = 0    # DATCOM assumption
    Y_r_wing = 0
    Y_r = Y_r_fin + Y_r_body + Y_r_wing

    # Y_zeta
    Y_zeta_fin = S_fin/S_w*a_r
    Y_zeta_rudder = a_v*ratio_alpha_dcl*alpha_dcl*K_*K_b*S_fin/S_w_total                   # demander pour formule !!
    return Y_v, Y_p, Y_r, Y_zeta_fin

def L_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    a_w = CL_alpha_wing(M)
    a_v = CL_beta_tail_vertical(M)
    a_r = CL_beta_rudder(M)

    # L_v
    grad_sidewash = sidewash_beta()
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_F = x_ac_v - x_cg
    V_F = (S_fin*l_F)/(S_w_total*b_w)
    beta = np.sqrt(1-M**2)
    
    Clbeta_CL_half_sweep = 0.0002
    KM_sweep = 1
    Clbeta_CL_A = 0
    Clbeta_dihedral = -0.00035
    KM_dihedral = 1.005
    L_v_fin = -V_F*z_F/l_F*a_v*(1-grad_sidewash)
    L_v_wing = 180/np.pi*( a_w*(Clbeta_CL_half_sweep*KM_sweep + Clbeta_CL_A) + dihedral_w*Clbeta_dihedral)
    L_v_body = (180*np.sqrt(AR_w)*d)/(np.pi*b_w)*(-0.0005*d/b_w-2.4*np.pi*Z_w/(180*b_w))
    L_v = L_v_fin + L_v_wing + L_v_body

    # L_p
    beta_Clp_k = -1.1       # approximatif !!
    k = a0_w/(2*np.pi)
    Clp_dihedral = 1 - 2*Z_w/(b_w/2)*np.sin(dihedral_w)+3*(Z_w/(b_w/2))**2*(np.sin(dihedral_w))**2  # C_lp_Gamma/C_lp_Gamma=0
    L_p_fin = -2*S_fin/S_w*a_v*(1-grad_sidewash)*z_F**2/b_v
    L_p_wb = beta_Clp_k*k/beta*Clp_dihedral 
    L_p = L_p_fin + L_p_wb

    # L_r
    B = np.sqrt(1-(M*np.cos(sweep_w))**2)
    L_r_fin = 2*S_fin/S_w*a_v*(1-grad_sidewash)*((z_F*np.sin(alpha)+l_F*np.cos(alpha))*(z_F*np.cos(alpha_e)-l_F*np.sin(alpha_e)))/b_v**2    
    L_r_body = 0    # DATCOM assumption
    clr_CL_num = 1+(AR_w*(1-B**2))/(2*B*(AR_w*B+2*np.cos(sweep_w))) + (AR_w*B+2*np.cos(sweep_w))/(AR_w*B+4*np.cos(sweep_w))*(np.tan(sweep_w)**2)/8
    clr_CL_den = 1+(AR_w+2*np.cos(sweep_w))/(AR_w+4*np.cos(sweep_w))*(np.tan(sweep_w)**2)/8
    clr_CL_M0 = 0.24
    clr_CL = (clr_CL_num/clr_CL_den)*clr_CL_M0
    dclr_dihedral = 1/12*(np.pi*AR_w*np.sin(sweep_w))/(AR_w+4*np.cos(sweep_w))
    L_r_wing = a_w*clr_CL + dclr_dihedral*dihedral_w
    L_r = L_r_fin + L_r_body + L_r_wing

    # L_zeta
    L_zeta_fin = V_F*z_F/l_F*a_r

    return L_v, L_p, L_r, L_zeta_fin

def N_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    beta = np.sqrt(1-M**2)
    a_w = CL_alpha_wing(M)
    a_v = CL_beta_tail_vertical(M)
    a_r = CL_beta_rudder(M)
    grad_sidewash = sidewash_beta()

    # N_v
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_F = x_ac_v - x_cg
    V_F = (S_fin*l_F)/(S_w_total*b_w)
    N_v_fin = V_F*a_v*(1-grad_sidewash)
    N_v_wing = 0 # neglected by DATCOM
    S_Bs = 0.51682
    K_N = 0.0042
    Re_lfus = V0*l_fus*rho/nu
    K_Rl = 1.52
    N_v_body = -180/np.pi*(K_N*K_Rl*S_Bs*l_fus)/(S_w_total*b_w)
    N_v = N_v_fin + N_v_wing + N_v_body

    # N_p
    beta_Clp_k = -1.1       # approximatif !!
    k = a0_w/(2*np.pi)
    Clp_dihedral = 1 - 2*Z_w/(b_w/2)*np.sin(dihedral_w)+3*(Z_w/(b_w/2))**2*(np.sin(dihedral_w))**2  # C_lp_Gamma/C_lp_Gamma=0
    L_p_fin = -2*S_fin/S_w*a_v*(1-grad_sidewash)*z_F**2/b_v
    L_p_wb = beta_Clp_k*k/beta*Clp_dihedral 
    L_p = L_p_fin + L_p_wb
    B = np.sqrt(1-M**2*(np.cos(sweep_w))**2)

    h0 = np.abs(x_ac_w-x_debut_wing)/c_mac_w
    h = np.abs(x_cg-x_debut_wing)/c_mac_w
    cnp_CL_M0 = -1/6*(AR_w+6*(AR_w+np.cos(sweep_w)))/(AR_w+4*np.cos(sweep_w))*((h0-h)*np.tan(sweep_w)/AR_w+((np.tan(sweep_w))**2)/12)
    cnp_CL = (AR_w+4*np.cos(sweep_w))/(AR_w*B+4*np.cos(sweep_w))*(AR_w*B+0.5*(AR_w*B+np.cos(sweep_w))*(np.tan(sweep_w))**2)/(AR_w+0.5*(AR_w+np.cos(sweep_w))*(np.tan(sweep_w))**2)*cnp_CL_M0
    N_p_fin = 2*S_fin/S_w*a_v*(1-grad_sidewash)*((z_F*np.sin(alpha)+l_F*np.cos(alpha))*(z_F*np.cos(alpha_e)-l_F*np.sin(alpha_e)))/b_v**2  
    N_p_body = 0 # assumption of DATCOM
    N_p_wing = -L_p_wb*np.tan(alpha)+L_p*np.tan(alpha)+cnp_CL*a_w
    N_p = N_p_fin + N_p_body + N_p_wing

    # N_r
    N_r_fin = -2*S_fin/S_w*a_v*(1-grad_sidewash)*(z_F*np.sin(alpha)+l_F*np.cos(alpha))**2/b_v**2
    N_r_body = 0
    cnr_CL2 = -0.03
    cnr_CD0 = -0.3
    N_r_wing = cnr_CL2*a_w**2 + cnr_CD0*CD0
    N_r = N_r_fin + N_r_body + N_r_wing

    # N_zeta
    #N_zeta_fin = -V_F*a_r*zeta ?

    return N_v, N_p, N_r


N_derivatives(W_b,x_b,W_crew1,W_crew2,V0)