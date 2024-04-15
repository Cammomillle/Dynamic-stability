import numpy as np 
from matplotlib import pyplot as plt 
import control as ct 
import pos_ellipse as pe
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit, fmin
from data import *
from functions import *

"""
THE FOLLOWING VALUES MUST CHANGE !!
"""
CD0=0.00964
T=-4.8+273.15 # temperature at 10 000 ft ! ->    shit ass guess !!!!
nu = 1.694*10**(-5) #Ns/m^2
rho = 0.909 # at 10,000 feet (assumed equilibrium position)
V0 = 28.809 # m/s
gamma = 1.4 

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
    
"------------------------------------------------Matrices----------------------------------------"
def compute_lat_matrices(m,Y_p,Y_r,Y_v,Ue,We,theta_e,L_p,L_r,L_v,N_p,N_r,N_v):
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
    
    """C1=np.matrix([[Y_ksi,Y_zeta],
                 [L_ksi,L_zeta],
                 [N_ksi,N_zeta],
                 [0,0],
                 [0,0]])"""

    A_1=-np.linalg.inv(A1) @ B1
    #B_1=np.linalg.inv(A1) @ C1
    return A_1

def compute_long_matrices(m,theta_e,Iy,X_w_dot,Z_w_dot,M_w_dot,X_u,X_w,X_q,Z_u,Z_w,Z_q,M_u,M_w,M_q,U_e,W_e):
    A2=np.matrix([[m, -X_w_dot, 0, 0],
                 [0, (m-Z_w_dot), 0, 0],
                 [0, -M_w_dot, Iy, 0],
                 [0, 0, 0, 1]])
    
    B2=np.matrix([[-X_u, -X_w, -(X_q-m*W_e), m*g*np.cos(theta_e)],
                 [-Z_u, -Z_w, -(Z_q+m*U_e), m*g*np.sin(theta_e)],
                 [-M_u, -M_w, -M_q, 0], 
                 [0, 0, -1, 0]])
    
    """C2=np.matrix([[X_eta, X_tau], 
                  [Z_eta, Z_tau],
                  [M_eta, M_tau], 
                  [0, 0]])"""

    A_2=-np.linalg.inv(A2) @ B2
    #B_2=np.linalg.inv(A2) @ C2
    return A_2

"-------------------------------------------------Utilities--------------------------------------"
def CL_alpha_wing(M):
    beta=np.sqrt(1-M**2)
    k=a0_w/(2*np.pi)
    CL_alphaW=(2*np.pi*AR_w)/(2+np.sqrt((AR_w*beta/k)**2*(1+np.tan(sweep_w_half)**2/beta**2)+4))
    return CL_alphaW

def CL_alpha_horizontal_tail(M):
    beta=np.sqrt(1-M**2)
    k=dclh_alpha_h/(2*np.pi)
    CL_alphaH=(2*np.pi*AR_h)/(2+np.sqrt((AR_h*beta/k)**2*(1+np.tan(sweep_h_half)**2/beta**2)+4))
    return CL_alphaH

def CL_alpha_vertical_tail(M):
    beta=np.sqrt(1-M**2)
    k=dclF_dbeta/(2*np.pi)
    AR_vBv_ratio = 1.02
    AR_vHBvB_ratio = 1.7
    Kh = 1
    AR_v_eff = (AR_vBv_ratio)*AR_v*(1+Kh*(AR_vHBvB_ratio-1))
    CL_alphaV=(2*np.pi*AR_v_eff)/(2+np.sqrt((AR_v_eff*beta/k)**2*(1+np.tan(sweep_v_half)**2/beta**2)+4)) 
    return CL_alphaV

def compute_CL(V_0):
    return 2*W/(rho*V_0**2*S_w_total)

def compute_CL_bis(M, alpha):
    K_WB = 1 - 0.25*(d_fus/b_w)**2 + 0.025*(d_fus/b_w)
    CL_alpha_WB = K_WB*CL_alpha_wing(M)
    grad_downwash = downwash(M)
    CL_alpha = CL_alpha_WB + CL_alpha_horizontal_tail(M)*eta_h*S_h/S_w_total*(1-grad_downwash)                                       
    CL = CL_alpha*(alpha-alpha_0)
    return CL
    
def compute_response(A,B,C,D,input_,T,X0):
    sys=ct.StateSpace(A,B,C,D)
    response=ct.forced_response(sys,T,input_,X0)

"----------------------------------------------DATCOM longitudinal derivatives--------------------------------"
def alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    " SLIDE 19 DATCOM "
    a=np.sqrt(R*gamma*T)
    Mach=V0/a 
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)

    #---- CL derivative ----
    CL=compute_CL(V0)
    K_WB=1-0.25*(d_fus/b_w)**2+0.025*(d_fus/b_w)
    grad_downwash = downwash(Mach)
    CL_alpha_WB=K_WB*CL_alpha_wing(Mach) #KWB*Cl_alpha_wing
    CL_alpha=CL_alpha_WB+CL_alpha_horizontal_tail(Mach)*eta_h*S_h/S_w_total*(1-grad_downwash)

    #--- CM derivative ----
    M_alpha = 0.5*rho*V0**2 / 36.5 * body_sum(Mach, x_cg)
    d_xac = -M_alpha / (0.5*rho*V0**2 * S_w_total * c_mac_w * CL_alpha_wing(Mach))
    x_acwb = x_ac_w/c_mac_w + d_xac
    num = x_acwb + CL_alpha_horizontal_tail(Mach)/CL_alpha_WB * eta_h * S_h/S_w_total * x_ac_h * (1-grad_downwash)  # diviser x_ac par c_mac ? 
    den = 1 + CL_alpha_horizontal_tail(Mach)/CL_alpha_WB * eta_h * S_h/S_w_total * (1-grad_downwash)
    x_ac = num/den
    Kn = -(x_cg/c_mac_w-x_ac)
    #print(x_ac, x_ac_w)
    #print('Stability margin', K)
    CM_alpha=-Kn*CL_alpha

    #--- CD derivative ----
    CD_alpha=2*CL*CL_alpha/(np.pi*AR_w*e_w)

    return CD_alpha, CL_alpha, CM_alpha


def u_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha,alpha_0):
    a=np.sqrt(R*gamma*T)
    Mach=V0/a 

    #---- CL derivative ----
    CL_alpha_for_mach,CD_alpha,CM_alpha=alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0) #CL_alpha for a given mach M !
    #CL_for_mach=CL_alpha_for_mach*(alpha-alpha_0)
    #CL_u=Mach**2/(1-Mach**2) * CL_for_mach # ATTENTION recompute CL for the given mach !!

    # Calcul alternatif donné par Prof. Dimitriadis
    CL = compute_CL_bis(Mach, alpha)                     # alpha = alpha_e assumption !!!
    CL_u = 2*CL - alpha_e*CL_alpha_for_mach
    
    #---- CD derivative ----
    #V0_plus=V0+2 #an increment on the speed
    #Mach_plus=V0_plus/a
    #dM=Mach_plus-Mach
    #CL_alpha_plus,_,_=alpha_derivatives(W_b, x_b, W_crew1, W_crew2, V0_plus)
    #CL_plus=CL_alpha_plus*(alpha-alpha_0)
    #dCD_dM=(CL_plus**2-CL_for_mach**2)/(dM*np.pi*AR_w*e_w) #We assumed that CD0 was independent of M !! -> is it okay ?
    "dCD_dm doit être calculé comme suit : 1) calculer CL pour M et pour M+delta M et en tirer les CD associés depuis la drag polar puis déf de la dérivée "
    #CD_u=Mach*dCD_dM
    "Here we will assume that x_ac_w doesn't move with the mach, which should be fine for low Mach n# (don't know any better)"

    # Calcul alternatif donné par Prof. Dimitriadis
    CD = CD0 + CL**2/(np.pi*AR_w*e_w)
    CD_u = 2*CD - alpha_e*CD_alpha

    #---- CM derivative ----
    #x_ac_plus=x_w
    #x_ac_for_mach=x_w 
    #dxac_w_dM=x_ac_plus-x_ac_for_mach/dM
    #CM_u=-CL_for_mach * dxac_w_dM

    # Calcul alternatif donné par Prof. Dimitriadis
    CM = 0  # at equilibrium                                       
    CM_u = 2*CM-alpha_e*CM_alpha

    return CD_u, CL_u, CM_u

def q_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    a=np.sqrt(R*gamma*T)
    M=V0/a 
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    Xw=(x_cg - x_ac_w)   # "taken positive rearward"
    
    #---- CD derivative ----
    CD_q=0 # DATCOM recommends to neglect it
    
    #---- CL derivative ----
    B = np.sqrt(1-(M*np.cos(sweep_w))**2)
    CL_qW = (AR_w+2*np.cos(sweep_w))/(AR_w*B+2*np.cos(sweep_w))*(1/2+2*Xw/c_mac_w)*CL_alpha_wing(M=0)
    l_T = x_ac_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)
    CL_qH = 2*CL_alpha_horizontal_tail(M)*eta_h*V_t_bar
    CL_q = CL_qW+CL_qH
    
    #---- CM derivative ----
    num = (AR_w**3*(np.tan(sweep_w))**2)/(AR_w*B+6*np.cos(sweep_w)) + 3/B   # B est une correction de compressibilité => B ~ 1
    den = (AR_w**3*(np.tan(sweep_w))**2)/(AR_w*B+6*np.cos(sweep_w)) + 3
    num1 = AR_w*(2*((Xw/c_mac_w)**2)+0.5*(Xw/c_mac_w))/(AR_w+2*np.cos(sweep_w))
    num2 = 1/24*AR_w**3*(np.tan(sweep_w))**2/(AR_w+6*np.cos(sweep_w)) 
    num3 = 1/8
    
    K = 0.9 # voir graphe slide 42 DATCOM
    CM_qw0 = - K*a0_w*np.cos(sweep_w)*(num1 + num2 + num3)
    CM_qw = CM_qw0*(num/den)
    CM_qh = -2*CL_alpha_horizontal_tail(M)*eta_h*V_t_bar*l_T/c_mac_w 
    print("\n Xw:",Xw, "CL_alpha_horizontal_tail", CL_alpha_wing(M), "V_t_bar", V_t_bar, "l_T", l_T, "\n")

    CM_q = CM_qw + CM_qh
    print("CMqh:", CM_qh, "CMqw:", CM_qw)

    return CD_q, CL_q, CM_q

def alpha_dot_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    # all negligible as advised by Prof. Dimitriadis
    #---- CD derivative ----
    CD_alpha_dot = 0 

    #---- CL derivative ----
    CL_alpha_dot = 0

    #---- CM derivative ----
    CM_alpha_dot = 0

    return CD_alpha_dot, CL_alpha_dot, CM_alpha_dot

def downwash(M):
    K_A = 1/AR_w - 1/(1+AR_w**1.7)
    K_lambda = (10-3*lambda_w)/7
    K_H = (1-h_H/b_w)/((2*l_H/b_w)**(1/2))                                                  # h_H to check !!
    grad_down_M0 = 4.44*(K_A*K_lambda*K_H*np.cos(sweep_w))**(1.19)
    grad_down = grad_down_M0*CL_alpha_wing(M)/CL_alpha_wing(0)   
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
    
    global_dwash = downwash(M)
    l_H_bodysum = x_ac_h - (x_debut_wing+c_w_root)

    ret = 0.0
    for frame in frames:
        Wf = 2*pe.b_i(frame*1000)/1000

        if frame < x_debut_wing:
            xicf = (x_debut_wing - frame)/c_w_root
            dwash = hyperbole(xicf, *Rfigure11) * CL_alpha_wing(0) / (0.08 * 180 / np.pi)

        elif frame > x_debut_wing + c_w_root:
            dwash = (1 - global_dwash)*(frame - x_debut_wing - c_w_root)/l_H_bodysum
        else:
            dwash = 0.0

        ret += Wf**2 * dwash * dXi
    
    return ret

def long_derivatives(W_b,x_b,W_crew1,W_crew2,V0, U_e, W_e):
    CD_u, CL_u, CM_u = u_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha_e,alpha_0)
    CD_alpha, CL_alpha, CM_alpha = alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0)
    CD_q, CL_q, CM_q = q_derivatives(W_b,x_b,W_crew1,W_crew2,V0)
    CD_alpha_dot, CL_alpha_dot, CM_alpha_dot = alpha_dot_derivatives(W_b,x_b,W_crew1,W_crew2,V0)
    print("CD_u: ", CD_u, "CL_u: ", CL_u,"CM_u: ", CM_u)
    print("CD_alpha: ", CD_alpha, "CL_alpha: ", CL_alpha, "CM_alpha: ", CM_alpha)
    print("CD_q: ", CD_q, "CL_q: ", CL_q, "CM_q: ", CM_q)
    print("CD_alpha_dot: ", CD_alpha_dot, "CL_alpha_dot: ", CL_alpha_dot, "CM_alpha_dot: ", CM_alpha_dot)
    #CD_eta, CL_eta, CM_eta = eta_derivatives(V0)

    L_e = compute_total_mass(W_b, W_crew1, W_crew2)
    D_e = 170.378 #N
    C_ze = (-L_e*np.cos(alpha_e) - D_e*np.sin(alpha_e))/(0.5*rho*V0**2*S_w_total)
    C_xe = (L_e*np.sin(alpha_e) - D_e*np.cos(alpha_e))/(0.5*rho*V0**2*S_w_total)
    
    # must translate in dimensional values !
    X_u = CL_u*np.sin(alpha_e) - CD_u*np.cos(alpha_e)
    X_u = X_u*0.5*rho*V0*S_w_total  
    X_w = 1/(np.cos(alpha_e))*(-C_ze+CL_alpha*np.sin(alpha_e) - CD_alpha*np.cos(alpha_e))
    X_w = X_w*0.5*rho*V0*S_w_total          
    X_q = CL_q*np.sin(alpha_e) - CD_q*np.cos(alpha_e)
    X_q = X_q*0.5*rho*V0*S_w_total*c_mac_w    
    X_w_dot = 1/np.cos(alpha_e)*(CL_alpha_dot*np.sin(alpha_e)-CD_alpha_dot*np.cos(alpha_e))
    X_w_dot = X_w_dot*0.5*rho*S_w_total*c_mac_w

    Z_u = -(CL_u*np.cos(alpha_e) + CD_u*np.sin(alpha_e))
    Z_u = Z_u*0.5*rho*V0*S_w_total
    Z_w = 1/np.cos(alpha_e)*(C_xe - CL_alpha*np.cos(alpha_e) - CD_alpha*np.sin(alpha_e))
    Z_w = Z_w*0.5*rho*V0*S_w_total         
    Z_q = -(CL_q*np.cos(alpha_e) + CD_q*np.sin(alpha_e))
    Z_q = Z_q*0.5*rho*V0*S_w_total*c_mac_w
    Z_w_dot = -1/np.cos(alpha_e)*(CL_alpha_dot*np.cos(alpha_e) + CD_alpha_dot*np.sin(alpha_e))
    Z_w_dot = Z_w_dot*0.5*rho*S_w_total*c_mac_w

    M_u = CM_u
    M_u = M_u*0.5*rho*V0*S_w_total*c_mac_w
    M_w = 1/np.cos(alpha_e)*CM_alpha
    M_w = M_w*0.5*rho*V0*S_w_total*c_mac_w 
    M_q = CM_q
    M_q = M_q*0.5*rho*V0*S_w_total*c_mac_w**2
    M_w_dot = 1/np.cos(alpha_e)*CM_alpha_dot
    M_w_dot = M_w_dot*0.5*rho*S_w_total*c_mac_w**2

    #Z_eta = -CL_eta*np.cos(alpha_e) - CD_eta*np.sin(alpha_e)
    #X_eta = CL_eta*np.sin(alpha_e) - CD_eta*np.cos(alpha_e)
    #M_eta = CM_eta

    return X_u, X_w, X_q, X_w_dot, Z_u, Z_w, Z_q, Z_w_dot, M_u, M_q, M_w, M_w_dot

"----------------------------------------------DATCOM lateral derivatives--------------------------------"
def sidewash_beta():    # sidewash derivative wrt beta
    sigma_beta = -0.276 + 3.06*S_fin/S_w_total*1/(1+np.cos(sweep_v)) + 0.4*z_w/d + 0.009*AR_w # slide 13 L13
    return sigma_beta

def derivative(func):
    """Make a derivative function of a function passed as argument"""
    eps = 0.000001
    return lambda x: (func(x+eps) - func(x-eps))/(2*eps)

def compute_S0():
    # Assumption: the body is assumed to be a circular revolution of radius
    # defined by its lower line (pe.ll_i)

    # Find x1 (maximum of derivative of body shape; here the lower body line is used)
    dsdx = lambda x: -derivative(pe.ll_i)(x)
    x1 = fmin(dsdx, 3000)[0]/1000

    # Get x0 (Corr. slide 49)
    x0 = l_fus*(0.378 + 0.527*(x1/l_fus))

    # Get a and b from body shape at x0
    a = pe.a_i(x0*1000)/1000
    b = pe.b_i(x0*1000)/1000

    # Return ellipse area
    return np.pi * a * b

def compute_Ki():
    # Find d (maximum body height at wing-body intersection)
    def bh_at_wing_intersection(x):
        if x < x_debut_wing:
            return x_debut_wing - x
        elif x > x_debut_wing + c_w_root:
            return x - (x_debut_wing + c_w_root)
        else:
            return -pe.bh_i(x*1000)/1000
    xd = fmin(bh_at_wing_intersection, 3.3)[0]
    d_Ki = pe.bh_i(xd*1000)/1000
    r = z_w/(d_Ki/2)

    # Return K_i (Corr. extracted from slide 49)
    if r >= 0:
        # Low wing case
        return 1.4915*r
    else:
        # High wing case
        return -1.844*r

def compute_Kb(eta):
    # Get Kb (extracted from slide 65; for lambda = 0.5)
    Kb_lambda_5 = np.genfromtxt('LatData/Kb_lambda_5.csv', delimiter=',')
    Kb_lambda_5 = interp1d(Kb_lambda_5[:,0], Kb_lambda_5[:,1])
    return Kb_lambda_5(eta)

def compute_Kprime(df):
    # Get Kprime (extracted from slide 65; for cf/c = 0.28)
    Kprime_c_ratio_28 = np.genfromtxt('LatData/Kprime_c_ratio_28.csv', delimiter=',')
    Kprime_c_ratio_28 = interp1d(Kprime_c_ratio_28[:,0], Kprime_c_ratio_28[:,1])
    return Kprime_c_ratio_28(df)

def Y_betav(M):
    k = 1
    d_Ki = 0.924
    grad_sidewash_eta = 0.724 + 3.06*(S_fin/S_w_total)/(1+np.cos(sweep_v)) + 0.4*z_w/d_Ki + 0.009*AR_w
    return -k*CL_alpha_vertical_tail(M)*grad_sidewash_eta*S_fin/S_w_total

def Y_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    a_v = CL_alpha_vertical_tail(M)
    grad_sidewash = sidewash_beta()
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_V = x_ac_v - x_cg   

    # Y_v
    Y_v_fin = Y_betav(M)      
    Y_v_wing = -0.0001*np.abs(dihedral_w)*180/np.pi                         # Valeur modifiée avec ajout de *180/np.pi  
    S0 = compute_S0()
    Ki = compute_Ki()
    Y_v_body = -2*Ki*S0/S_w_total
    Y_v = Y_v_fin + Y_v_body + Y_v_wing
    print("Y_v: ", Y_v)
    Y_v = Y_v*0.5*rho*V0*S_w_total # dimensional

    # Y_p (roll rate)           
    Y_p_fin = -2*S_fin/S_w_total*a_v*(1-grad_sidewash)*(z_V*np.cos(alpha_e)-l_V*np.sin(alpha_e))/b_w
    Y_p_wing = 0 # assumed negligible
    Y_p_body = 0
    Y_p = Y_p_fin + Y_p_wing + Y_p_body
    print('Y_p: ', Y_p)
    Y_p = Y_p*0.5*rho*V0*S_w_total

    # Y_r (yaw rate)
    Y_r_fin = -2/b_w*(l_V*np.cos(alpha)+z_V*np.sin(alpha))*Y_betav(M)
    Y_r_body = 0    # DATCOM assumption
    Y_r_wing = 0
    Y_r = Y_r_fin + Y_r_body + Y_r_wing
    print("Y_r: ", Y_r)
    Y_r = Y_r*0.5*rho*V0*S_w_total

    # Y_zeta (rudder deflection)
    """eta = 20*np.pi/180                                   # eta = max elevator deflection ?
    df = 30*np.pi/180                                       # df = max rudder deflection ? pq flaps deflection dans la formule ?
    Kb = compute_Kb(eta)
    Kprime = compute_Kprime(df)
    Y_zeta = a_v*alpha_dcl_ratio*alpha_dcl*Kprime*Kb*S_fin/S_w_total"""          

    # Y_ksi (aileron deflection)
    """Y_ksi = 0   # negligible according to DATCOM"""

    return Y_v, Y_p, Y_r

def L_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    a_w = CL_alpha_wing(M)
    a_v = CL_alpha_vertical_tail(M)

    # L_v - pas modifié
    grad_sidewash = sidewash_beta()
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_V = x_ac_v - x_cg
    V_F = (S_fin*l_V)/(S_w_total*b_w)
    beta = np.sqrt(1-M**2)
    
    Clbeta_CL_half_sweep = 0.0002               
    KM_sweep = 1
    Clbeta_CL_A = 0
    Clbeta_dihedral = -0.00035
    KM_dihedral = 1.005
    L_v_fin = -V_F*z_V/l_V*a_v*(1-grad_sidewash)
    L_v_wing = 180/np.pi*( a_w*(Clbeta_CL_half_sweep*KM_sweep + Clbeta_CL_A) + dihedral_w*Clbeta_dihedral)
    L_v_body = (180*np.sqrt(AR_w)*d)/(np.pi*b_w)*(-0.0005*d/b_w-2.4*np.pi*z_w/(180*b_w))
    L_v = L_v_fin + L_v_wing + L_v_body
    print("L_v: ", L_v)
    L_v = L_v*0.5*rho*V0*S_w_total*b_w

    # L_p
    beta_Clp_k = -1.1       # approximatif !!
    k = a0_w/(2*np.pi)
    Clp_dihedral = 1 - 2*z_w/(b_w/2)*np.sin(dihedral_w)+3*(z_w/(b_w/2))**2*(np.sin(dihedral_w))**2  # C_lp_Gamma/C_lp_Gamma=0
    L_p_fin = -2*S_fin/S_w_total*a_v*(1-grad_sidewash)*z_V**2/b_w
    L_p_wb = beta_Clp_k*k/beta*Clp_dihedral 
    L_p = L_p_fin + L_p_wb
    print("L_p: ", L_p)
    L_p = L_p*0.5*rho*V0*S_w_total*b_w
    #L_p = L_p*0.25*rho*V0*S_w_total*b_w**2

    # L_r
    B = np.sqrt(1-(M*np.cos(sweep_w))**2)
    #L_r_fin = 2*S_fin/S_w_total*a_v*(1-grad_sidewash)*((z_V*np.sin(alpha)+l_F*np.cos(alpha))*(z_V*np.cos(alpha_e)-l_F*np.sin(alpha_e)))/b_w**2    
    L_r_fin = -2/b_w**2*(l_V*np.cos(alpha)+z_V*np.sin(alpha))*(z_V*np.cos(alpha)-l_V*np.sin(alpha))*Y_betav(M)

    L_r_body = 0    # DATCOM assumption

    clr_CL_num = 1+(AR_w*(1-B**2))/(2*B*(AR_w*B+2*np.cos(sweep_w))) + (AR_w*B+2*np.cos(sweep_w))/(AR_w*B+4*np.cos(sweep_w))*(np.tan(sweep_w)**2)/8
    clr_CL_den = 1+(AR_w+2*np.cos(sweep_w))/(AR_w+4*np.cos(sweep_w))*(np.tan(sweep_w)**2)/8
    clr_CL_M0 = 0.25
    clr_CL = (clr_CL_num/clr_CL_den)*clr_CL_M0
    dclr_dihedral = 1/12*(np.pi*AR_w*np.sin(sweep_w))/(AR_w+4*np.cos(sweep_w))
    C_L = compute_CL_bis(M, alpha)
    L_r_wing = C_L*clr_CL + dclr_dihedral*dihedral_w
    L_r = L_r_fin + L_r_body + L_r_wing
    print("L_r: ", L_r)
    L_r = L_r*0.5*rho*V0*S_w_total*b_w
    #L_r = L_r*0.25*rho*V0*S_w_total*b_w**2

    # L_zeta
    """Kb = compute_Kb(eta)
    Kprime = compute_Kprime(df)
    L_zeta = a_v*alpha_dcl_ratio*alpha_dcl*Kprime*Kb*S_fin/S_w_total*(z_V*np.cos(alpha_e)-l_F*np.sin(alpha_e))/b_w"""
    
    # L_ksi 
    """cld_theory = 3.7
    Kprime = 0.86
    cld_ratio = 
    # Linear interpolation to find betacl_k 
    betacl_k = (((0.8-0.44)/(8-4))*AR_w+0.44) - (((0.28-0.085)/(8-4))*AR_w+0.085)

    clalpha = 6.31206

    cldR = np.abs(alphadR)*cldprime_R
    cldL = np.abs(alphadL)*cldprime_L
    LdA = cldL+cldR"""

    return L_v, L_p, L_r

def N_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha):
    a = np.sqrt(R*gamma*T)
    M = V0/a
    beta = np.sqrt(1-M**2)
    a_w = CL_alpha_wing(M)
    a_v = CL_alpha_vertical_tail(M)
    grad_sidewash = sidewash_beta()
    x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    l_V = x_ac_v - x_cg

    # N_v - pas modifié
    V_F = (S_fin*l_V)/(S_w_total*b_w)
    N_v_fin = V_F*a_v*(1-grad_sidewash)
    N_v_wing = 0 # neglected by DATCOM
    S_Bs = 0.51682
    K_N = 0.0042
    Re_lfus = V0*l_fus*rho/nu
    K_Rl = 1.52
    N_v_body = -180/np.pi*(K_N*K_Rl*S_Bs*l_fus)/(S_w_total*b_w)
    N_v = N_v_fin + N_v_wing + N_v_body
    print("N_v: ", N_v)
    N_v = N_v*0.5*rho*V0*S_w_total*b_w

    # N_p
    beta_Clp_k = -1.1       # approximatif !!
    k = a0_w/(2*np.pi)
    Clp_dihedral = 1 - 2*z_w/(b_w/2)*np.sin(dihedral_w)+3*(z_w/(b_w/2))**2*(np.sin(dihedral_w))**2  # C_lp_Gamma/C_lp_Gamma=0
    L_p_fin = -2*S_fin/S_w_total*a_v*(1-grad_sidewash)*z_V**2/b_w
    L_p_wb = beta_Clp_k*k/beta*Clp_dihedral 
    L_p = L_p_fin + L_p_wb
    B = np.sqrt(1-M**2*(np.cos(sweep_w))**2)

    h0 = np.abs(x_ac_w-x_debut_wing)/c_mac_w
    h = np.abs(x_cg-x_debut_wing)/c_mac_w
    cnp_CL_M0 = -1/6*(AR_w+6*(AR_w+np.cos(sweep_w)))/(AR_w+4*np.cos(sweep_w))*((h0-h)*np.tan(sweep_w)/AR_w+((np.tan(sweep_w))**2)/12)
    cnp_CL = (AR_w+4*np.cos(sweep_w))/(AR_w*B+4*np.cos(sweep_w))*(AR_w*B+0.5*(AR_w*B+np.cos(sweep_w))*(np.tan(sweep_w))**2)/(AR_w+0.5*(AR_w+np.cos(sweep_w))*(np.tan(sweep_w))**2)*cnp_CL_M0
    N_p_fin = 2*S_fin/S_w_total*a_v*(1-grad_sidewash)*((z_V*np.sin(alpha)+l_V*np.cos(alpha))*(z_V*np.cos(alpha_e)-l_V*np.sin(alpha_e)))/b_w**2  
    N_p_body = 0 # assumption of DATCOM
    N_p_wing = -L_p_wb*np.tan(alpha)+L_p*np.tan(alpha)+cnp_CL*a_w
    N_p = N_p_fin + N_p_body + N_p_wing
    print("N_p: ", N_p)
    N_p = N_p*0.5*rho*V0*S_w_total*b_w
    #N_p = N_p*0.25*rho*V0*S_w_total*b_w**2

    # N_r
    N_r_fin = 2/b_w**2*(l_V*np.cos(alpha)+z_V*np.sin(alpha))**2*Y_betav(M)
    N_r_body = 0
    cnr_CL2 = -0.03
    cnr_CD0 = -0.3
    C_L = compute_CL_bis(M, alpha)
    N_r_wing = cnr_CL2*C_L**2 + cnr_CD0*CD0
    N_r = N_r_fin + N_r_body + N_r_wing
    print("N_r: ", N_r)
    N_r = N_r*0.5*rho*V0*S_w_total*b_w
    #N_r = N_r*0.25*rho*V0*S_w_total*b_w**2

    # N_zeta
    """Kb = compute_Kb(eta)
    Kprime = compute_Kprime(df)
    N_zeta = -a_v*alpha_dcl_ratio*alpha_dcl*Kprime*Kb*S_fin/S_w_total*(z_V*np.sin(alpha_e)+l_F*np.cos(alpha_e))/b_w"""    

    # N_ksi
    return N_v, N_p, N_r

def long_modes_caract(eigenVals): # gives the natural frequency and damping ratio of the longitudinal modes (phugoid and short period)
    # Imaginary part correspond to damped natural frequency : higher freq = short period oscillations and lower freq = phugoids
    phugoids_index = np.where(eigenVals.imag == min(abs(eigenVals.imag)))
    short_period_index = np.where(eigenVals.imag == max(abs(eigenVals.imag)))
    phugoids = complex(eigenVals[phugoids_index][0])
    short_period = complex(eigenVals[short_period_index][0])

    # Natural frequencies computation
    omega_n_phug = np.sqrt((phugoids.real)**2 + (phugoids.imag)**2)         # phugoids 
    omega_n_sp = np.sqrt((short_period.real)**2 + (short_period.imag)**2)   # short period 

    # Damping ratios computation
    damp_phug = -phugoids.real/omega_n_phug
    damp_sp = - short_period.real/omega_n_sp

    print("omega_n_phug:", omega_n_phug, "damp_phug:", damp_phug, "omega_n_sp:", omega_n_sp, "damp_sp:", damp_sp)
    
    return omega_n_phug, damp_phug, omega_n_sp, damp_sp

def lat_modes_caract(eigenVals): # gives the natural frequency and damping ratio of the lateral modes (dutch roll, roll subsidence, spiral mode)
    # Complex conjudated values correspond to dutch roll, the most negative value correspond to roll subsidence and the lower amplitude one to the spiral mode
    dutch_roll_index = np.where(eigenVals.imag < 0) # only complex conjugated pair
    dutch_roll = complex(eigenVals[dutch_roll_index][0])
    roll_sub_index = np.where((eigenVals.real == min(eigenVals.real)) & (eigenVals.imag == 0))   # roll subsidence is highly damped
    roll_sub = float(eigenVals[roll_sub_index][0].real)
    spiral_index = np.where(eigenVals.real > 0)  # spiral mode is hardly damped
    spiral = float(eigenVals[spiral_index][0].real)

    # Natural frequencies and damping computation
    omega_n_dutch = np.sqrt((dutch_roll.real)**2 + (dutch_roll.imag)**2)
    damp_dutch = -dutch_roll.real/omega_n_dutch

    # Time constant computation (is formula tau = 1/(damp*omega_n) ok?)
    tau_roll_sub = 1/(np.abs(roll_sub.real))
    tau_spiral = 1/(np.abs(spiral.real))

    print("omega_n_dutch:", omega_n_dutch, "damp_dutch:", damp_dutch, "tau_roll_sub:", tau_roll_sub, "tau_spiral:", tau_spiral)

    return omega_n_dutch, damp_dutch, tau_roll_sub, tau_spiral

#-------- Longitudinal matrix A computation ---------
W_tot = compute_total_mass(W_b, W_crew1, W_crew2)
m = W_tot/g
theta_e = 0
U_e = V0*np.cos(alpha_e)
W_e = V0*np.sin(alpha_e)

X_u, X_w, X_q, X_w_dot, Z_u, Z_w, Z_q, Z_w_dot, M_u, M_q, M_w, M_w_dot = long_derivatives(W_b,x_b,W_crew1,W_crew2,V0, U_e, W_e)
A_long = compute_long_matrices(m,theta_e,Iy,X_w_dot,Z_w_dot,M_w_dot,X_u,X_w,X_q,Z_u,Z_w,Z_q,M_u,M_w,M_q,U_e,W_e)
eigenval_A_long, eigenvect_A_long = np.linalg.eig(A_long)
#print("A_long:", A_long)
#print("eigenVect of A_long: ", eigenvect_A_long)
print("eigenVal of A_long: ", eigenval_A_long)

long_modes_caract(np.array(eigenval_A_long))

#-------- Lateral matrix A computation --------
alpha = alpha_e  # aircraft aoa, ok ?

Y_v, Y_p, Y_r = Y_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha)
L_v, L_p, L_r = L_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha)
N_v, N_p, N_r = N_derivatives(W_b,x_b,W_crew1,W_crew2,V0, alpha)
A_lat = compute_lat_matrices(m,Y_p,Y_r,Y_v,U_e,W_e,theta_e,L_p,L_r,L_v,N_p,N_r,N_v)
eigenval_A_lat, eigenvect_A_lat = np.linalg.eig(A_lat)
#print("\n A_lat:", A_lat)
#print("eigenVect of A_lat: ", eigenvect_A_lat)
print("eigenVal of A_lat: ", eigenval_A_lat)

lat_modes_caract(np.array(eigenval_A_lat))

