import numpy as np 
from matplotlib import pyplot as plt 
import control as ct 
from data import *
from functions import *
import aerodynamic

"Coucou izlou"
"""
THE FOLLOWING VALUES ARE VERY SUSPECT THEY MUST DEFINITELY CHANGE !!
"""
CD0=0.014


"""  END OF SUSPECT VALUES """
g=9.81
W_b_w = 0*g     # ballasts at wing
W_b_t = 0*g     # ballasts at tail
W_b = W_b_w + W_b_t  # ballasts total weight

W_crew1 = 80*g   # 1 st crew
W_crew2 = 80*g   # 2 nd crew

#******* CG of the ballasts **********
x_b = 0
if(W_b!=0):
    x_b = (W_b_w*x_b_w + W_b_t*x_b_t)/W_b
T=-4.8+273.15 #temperature at 10000 ft ! -> shit ass guess !
R=287 #J/kg/k pour l'air
W=compute_total_mass(W_b,W_crew1,W_crew2)
2*W/(rho*V0**2*S_w_total)
print("CL",CL)
print("CL de aho",)
"""
   Comments: what about alpha, alpha_0 in u_derivatives ?  
             Add the swept angle lambda at c/4 and c/2 somewhere in the beginning of the code !
             h-h0 in q_derivatives must still be computed
"""
"-------------------------------------------------Utilities--------------------------------------"

def compute_CL(V_0):
    return 2*W/(rho*V_0**2*S_w_total)
def compute_x_ac(V_0):
    return coucou
def compute_response(A,B,C,D,input_,T,X0):
    sys=ct.StateSpace(A,B,C,D)
    response=ct.forced_response(sys,T,input_,X0)
    
"------------------------------------------------Matrices----------------------------------------"

def compute_matrices(m,Ix,Ixz,Iz,Y_p,Y_r,Y_v,Ue,We,theta_e,L_p,L_r,L_v,N_p,N_r,N_v,Y_ksi,Y_zeta,L_ksi,L_zeta,N_ksi,N_zeta):
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
    
    A=-np.linalg.inv(A1) @ B1
    B=np.linalg.inv(A1) @ C1
    
    
"----------------------------------------------Inertias------------------------------------------"    

def compute_inertias():
    return coucou

"----------------------------------------------DATCOM derivatives--------------------------------"
def CL_alpha_wing(M):
    beta=np.sqrt(1-M**2)
    k=a_w/2*np.pi 
    CL_alphaW=(2*np.pi*AR_W)/(2+np.sqrt((AR_W*beta/k)**2*(1+np.tan(Lambda_c_half)**2/beta**2)+1))
    return CL_alphaW

def CL_alpha_tail_horizontal(M):
    beta=np.sqrt(1-M**2)
    k=dclh_alpha_h/2*np.pi 
    CL_alphaH=(2*np.pi*AR_V)/(2+np.sqrt((AR_V*beta/k)**2*(1+np.tan(Lambda_c_half)**2/beta**2)+1))
    return CL_alphaH
def alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0):
    " SLIDE 19 DATCOM "
    " DIAMETRE FUSELAGE PUE DU CUL ON A PRIS EPAISSEUR FUSELAGE A LA PLACE MAIS CEST NUL!!!!"
    
    #CL derivative
    CL=compute_CL(V0)
    m = 2*(z_tail-z_wing)/b_w
    grad_downwash = 1.75*a_w/(np.pi * AR_w * ((2*lambda_w*l_T)/b_w)**0.25*(1+abs(m)))
    d=width_fus
    K_WB=1-0.25*(d/b_w)**2+0.025*(d/b_w)
    CL_alpha_WB=K_WB*a_w #KWB*Cl_alpha_wing 
    eta_H=0.95 # between 0.9 and 1 in slide 19 
    CL_alpha=CL_alpha_WB+dCLh_alpha_h*eta_H*S_h/S_w_total*(1-grad_downwash)
    x_acwb=x_w
    x_ach=x_t_h
    num=x_acwb+dCLh_alpha_h/CL_alpha_WB * eta_H * S_h/S_w_total * x_ach *(1-grad_downwash)
    den=1+dCLh_alpha_h/CL_alpha_WB * eta_H * S_h/S_w_total * (1-grad_downwash)
    x_ac=num/den
    "On n'a pas tenu compte des effets du body -> check slide 22-23 DATCOM"
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    
    #CM derivative
    CM_alpha=(x_cg-x_ac)*CL_alpha
    #CD derivative
    CD_alpha=2*CL*CL_alpha/(np.pi*AR_w*e)
    return CL_alpha,CD_alpha,CM_alpha


def u_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha,alpha_0):
    a=np.sqrt(R*gamma*T)
    Mach=V0/a #mach number
    #CL derivative
    CL_alpha_for_mach,CD_alpha,CM_alpha=alpha_derivatives(W_b,x_b,W_crew1,W_crew2,V0) #CL_alpha for a given mach M !
       
    CL_for_mach=CL_alpha_for_mach*(alpha-alpha_0)
    CL_u=Mach**2/(1-Mach**2) * CL_for_mach # ATTENTION recompute CL for the given mach !!
    
    V0_plus=V0+2 #an increment on the speed
    Mach_plus=V0_plus/a
    dM=Mach_plus-Mach
    CL_alpha_plus,_,_=alpha_derivatives(W_b, x_b, W_crew1, W_crew2, V0_plus)
    CL_plus=CL_alpha_plus*(alpha-apha_0)
    dCD_dM=(CL_plus**2-CL_for_mach**2)/(dM*np.pi*AR_W*e) #We assumed that CD0 was independent of M !! -> is it okay ?
    "dCD_dm doit être calculé comme suit : 1) calculer CL pour M et pour M+delta M et en tirer les CD associés depuis la drag polar puis déf de la dérivée "
    #CD derivative
    CD_u=Mach*dCD_dM
    
    "Here we will assume that x_ac_w doesn't move with the mach, which should be fine for low Machs (don't know any better)"
    #CM derivative
    x_ac_plus=x_w
    x_ac_for_mach=x_w 
    dxac_w_dM=x_ac_plus-x_ac_for_mach/dM
    CM_u=-CL_for_mach * dxac_w_dM
    return CL_u,CD_u,CM_u
def q_derivatives(W_b,x_b,W_crew1,W_crew2,V0,alpha,alpha_0):
    a=np.sqrt(R*gamma*T)
    Mach=V0/a #mach number
    
    #CD derivative
    CD_q=0 #DATCOM recommends to neglect it
    
    #CL derivative
    "Lambda is the c/4 sweep angle -> ?"
    "(h-h0)c_mean_mean is the distance between the center of gravity and the wing AC"
    Lambda=lambda_w
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w) # x-location of the LE of the projected c_mac of the wing
    x_ac_w=x_w
    h_0 = (x_ac_w-x_le_wing)/c_mac_w
    x_cg=compute_x_cg(W_b,x_b,W_crew1,W_crew2)
    h=compute_h(x_cg)
    B=np.sqrt(1-M**2*cos(Lambda))
    CL_qW=(AR_W+2*cos(Lambda))/(AR_W*B+2*cos(Lambda))*(1/2+2*(h-h0))*CL_alpha_wing(M=0)
    eta_H=0.95 # between 0.9 and 1 in slide 19 
    l_T = x_ac_t_h - x_cg
    V_t_bar = (l_T*S_h)/(S_w_total*c_mac_w)
    CL_qH=2*CL_alpha_tail_horizontal(M=Mach)*eta_H*V_T_bar
    CL_q=CL_qW+CL_qH
    
    #CM derivative
    return coucou
    
