import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
from data import *
from functions import *
from matplotlib.ticker import FuncFormatter, AutoLocator

plt.rcdefaults()  # restore the defaults
from cycler import cycler

custom_colors = ['#00707f', '#f07f3c', '#7db928', '#ffd000', '#e62d31', '#5b57a2']
custom_cycler = cycler(color=custom_colors)
plt.rcParams.update({
    'axes.prop_cycle': custom_cycler,
    'text.usetex': True,
    'text.latex.preamble': '\\usepackage{stix2}',
    'font.family': 'serif',
    'font.serif': 'stix2'
})

#**********************************************************************************************
#*** Stability (pitch, yaw) for a configuration (number of pilot, pilots' weight, ballasts) *** 
#**********************************************************************************************

#******* Weights [kilograms * g] **********
W_b_crew = 0*g     # ballasts at wing
W_b_t = 0*g     # ballasts at tail
W_b = W_b_crew + W_b_t  # ballasts total weight

W_crew1 = 80*g   # 1 st crew
W_crew2 = 122*g   # 2 nd crew

#******* CG of the ballasts **********
x_b = 0
if(W_b!=0):
    x_b = (W_b_crew*x_b_w + W_b_t*x_b_t)/W_b
    
def round_formatter(value, pos):
    return round(value, 2) 

def plot_results(x_cg_enveloppe, h_n):

    x_n = h_n*c_mac_w
    x_cg_min = compute_x_cg(37.5*g, x_b_w, 40*g, 0) # CG position for 1 crew of 90lb
    x_cg_max = compute_x_cg(15*g, x_b_t, 122*g, 122*g) # CG position for 2 crews of 270lb each
    plt.scatter([x_cg_enveloppe[0]*3.28084,x_cg_enveloppe[1]*3.28084],[0,0],label="Static margin range",color="darkorange", s=80)
    plt.scatter(x_ac_w*3.28084,0,label="AC", s=80, color="yellowgreen")
    plt.scatter(x_n*3.28084,0,label="NP", s=80, color='red')
    plt.scatter([x_cg_min*3.28084, x_cg_max*3.28084], [0,0], label="CG variation", s=80, color="darkcyan")
    plt.plot([2.6*3.28084,3.7*3.28084],[0,0],color="black")
    """x_tab=np.linspace(-c__w/2,c__w/2,100)
    a=c__w/2
    b=0.5
    y_tab=[]
    x_cc=[]
    for x in x_tab :
        y_tab.append(np.sqrt(1-x**2/a**2)*b)
        x_cc.append(x)
    for x in x_tab:
        y_tab.append(-np.sqrt(1-x**2/a**2)*b)
        x_cc.append(x)
    x_cc=np.array(x_cc)
    plt.plot(x_cc+x_debut_wing,y_tab,label="wings")"""
    ax = plt.gca()
    ax.set_ylim(-0.5,0.5)
    ax.set_xlim((10.0,11.2))
    ticks=[x_cg_enveloppe[0],x_cg_enveloppe[1],x_ac_w, x_n, x_cg_min, x_cg_max]
    ticks=np.array(ticks)*3.28084
    ax.set_xticks(ticks)
    ax.get_yaxis().set_visible(False)
    ax.xaxis.set_major_formatter(FuncFormatter(round_formatter))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.xlabel("Longitudinal position [ft]", fontsize=11)
    plt.legend(loc='best',fontsize=11)
    plt.savefig("CGs.pdf",bbox_inches='tight')
    plt.show()

def x_cg_enveloppe(h_n): # computation of the CG enveloppe diagram
    """
    We want Kn to be such as : Kn>0.05 & Kn<0.25 or Kn=hn-h iif h=hn-K which means that h should be such as h>hn-0.25 and h<hn-0.05
    On peut prendre une plus grande marge
    """
    h_low = h_n-K_n_high  # h_min to be stable
    h_high = h_n-K_n_low  # h_max to be stable
    x_cg_low = c_mac_w*h_low
    x_cg_high = c_mac_w*h_high
    print("Safe CG enveloppe: ", "[", x_cg_low, ";", x_cg_high, "]", "[m] \n")
    return x_cg_low, x_cg_high
    
def ballast_positions():
            # On calcule le x_cg sans ballasts et on vérifie la stabilité
            W_b=0
            x_b=0
            #x_cg=compute_x_cg(W_b, x_b, W_crew1, W_crew2)
            #h=compute_h(x_cg)
            #h_n=compute_hn(x_cg)
            #K_n,is_kn_ok=compute_K_n(h_n,h)
            K_n,is_kn_ok=compute_K_n_bis(W_b, x_b, W_crew1, W_crew2)

            if is_kn_ok==True:
                print("Stable without ballasts")
                return x_cg
            if is_kn_ok==False:
              m_lim_tot=60 #Total ballast weight admissible
              m_lim_tail=15 #Maximum ballast weight at the tail 
              m_lim_crew=m_lim_tot-m_lim_tail #Maximum ballast weight under the crew
              mb1=np.arange(0,m_lim_crew+1,7.5)*g
              mb2=np.arange(0,m_lim_tail+1,7.5)*g
              fine=False
              best_list=[]
              W_b_list=[]
              for k in range(len(mb1)):
                  if(fine==True):
                      break
                  for l in range(len(mb2)):
                      W_b=mb1[k]+mb2[l]
                      if(W_b==0):
                          continue
                      x_b=(x_b_w*mb1[k]+x_b_t*mb2[l])/W_b # center of gravity of the ballasts 
                      x_cg_ballast=compute_x_cg(W_b, x_b, W_crew1, W_crew2)
                      #h=compute_h(x_cg_ballast)
                      #h_n=compute_hn(x_cg)
                      K_n,is_kn_ok=compute_K_n_bis(W_b, x_b, W_crew1, W_crew2)
                      if(is_kn_ok==True):
                          fine=True
                          W_b_list.append(W_b)
                          best_list.append([mb1[k],mb2[l],x_cg_ballast])
              try:          
                Wbmin=W_b_list[np.argmin(W_b_list)]
                print("[weight b1, weight b2, CG with ballasts]: ", best_list)
              except:
                print("No stable configurations with ballasts")
                return 0, 0
              for val in best_list:
                  print(val)
                  if val[0]+val[1]==Wbmin:
                      w_b_1=val[0]
                      w_b_2=val[1]
                      x_cg_ballast=val[2]
                      break
                          
              if(fine==False):
                  print("No stable configurations with ballasts")
                  return 0, 0
              else:
                  print("Ballasts: {0} kg under crew1 seat and ".format(w_b_1/g)+"{0} kg at the tail CG ".format(w_b_2/g))
                  return w_b_1+w_b_2, (w_b_1*x_b_w+w_b_2*x_b_t)/(w_b_1+w_b_2)

def yaw_stability(l_cg, l_F): # see Conceptual design slides 56-57
    dCn_i_beta = -0.017   # empirical values for high wings                  !!! à changer si on modifie la position des ailes !!!
    Kb = 0.3*l_cg/l_fus + 0.75*h_f_max/l_fus - 0.105
    dCn_fus_beta = -Kb*S_fuselage*l_fus/(S_w*b_w)*np.sqrt(hf1/hf2)*np.cbrt(bf2/bf1) # effect of the fuselage
    dCn_fin_beta = dCLF_dbeta * S_fin*l_F/(S_w*b_w) # effect of the fin
    dCn_beta = dCn_fus_beta + dCn_i_beta + dCn_fin_beta
    print("Fuselage:", dCn_fus_beta)
    print("Fin:", dCn_fin_beta)
    dCn_beta_is_ok = False
    if(dCn_beta>0):
        dCn_beta_is_ok = True # if dC_n/db > 0, yaw stability is good !
    return dCn_beta, dCn_beta_is_ok

# Total mass of the sailplane
total_mass = compute_total_mass(W_b, W_crew1, W_crew2)
#print("Total mass: ", total_mass/g, "kg")

# CG of the sailplane
x_cg = compute_x_cg(W_b, x_b, W_crew1, W_crew2)
#print("CG of the sailplane: ", x_cg, "m")
print("CG of the wing: ", x_w, "m")
print("CG of the  sailplane before ballasts computation: ", x_cg, "m")

# CG of the empty sailplane (i.e., without ballasts and crew)
x_cg_empty = compute_x_cg_empty()
print("CG of the empty sailplane: ", x_cg_empty, "m")

# Static pitch stability of the sailplane without ballasts
#h=compute_h(x_cg)
#h_n=compute_hn(x_cg)
#K_n=compute_K_n(h_n, h)
K_n, is_kn_ok=compute_K_n_bis(W_b, x_b, W_crew1, W_crew2)
h_n = K_n + x_cg/c_mac_w
print("Pitch stability Kn: ", K_n, "[-]")

# Enveloppe of the stable CG 
x_cg_env = x_cg_enveloppe(h_n)

if is_kn_ok == False:
    # CG with ballasts 
    W_ballasts, x_cg_ballasts = ballast_positions()

    # Static pitch stability of the sailplane with ballasts
    #h=compute_h(x_cg)
    #h_n=compute_hn(x_cg)
    #K_n=compute_K_n(h_n, h)
    K_n=compute_K_n_bis(W_ballasts, x_cg_ballasts, W_crew1, W_crew2)[0]
    h_n = K_n + x_cg/c_mac_w
    print("Pitch stability Kn: ", K_n, "[-]")

# Pitch stability diagram 
plot_results(x_cg_env,h_n)

# Yaw stability 
l_F = x_ac_v - x_cg # (see Conceptual Design slide 56)
yaw_stab = yaw_stability(x_cg, l_F)
print("Yaw stability dCn_beta:", yaw_stab, "[-]")


