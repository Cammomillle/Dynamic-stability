import numpy as np 
from matplotlib import pyplot as plt 
import control as ct 
from data import *
from functions import *
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
  
W=compute_total_mass(W_b,W_crew1,W_crew2)
CL=2*W/(rho*V_0**2*S_w_total) 
def compute_response(A,B,C,D,input_,T,X0):
    sys=ct.StateSpace(A,B,C,D)
    response=ct.forced_response(sys,T,input_,X0)

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
    
    
def compute_inertias():
    return coucou

def alpha_derivatives(W_b,x_b,W_crew1,W_crew2):
    " SLIDE 19 DATCOM "
    " DIAMETRE FUSELAGE PUE DU CUL ON A PRIS EPAISSEUR FUSELAGE A LA PLACE MAIS CEST NUL!!!!"
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
    C_M_alpha=(x_cg-x_ac)*CL_alpha
    C_D_alpha=2*CL*CL_alpha/(np.pi*AR_w*e)
def drag_derivatives():

def moment_derivativs():
    
    
