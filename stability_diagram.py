import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
from data import *
from functions import *
from cycler import cycler

"""plt.rcdefaults()  # restore the defaults

custom_colors = ['#00707f', '#f07f3c', '#7db928', '#ffd000', '#e62d31', '#5b57a2']
custom_cycler = cycler(color=custom_colors)
plt.rcParams.update({
    'axes.prop_cycle': custom_cycler,
    'text.usetex': True,
    'text.latex.preamble': '\\usepackage{stix2}',
    'font.family': 'serif',
    'font.serif': 'stix2'
})"""
#************************************************
#*** Pitch stability diagram of the sailplane *** 
#************************************************


def cg_diag():
    W_crew1=np.linspace(0,400,20) # lbs
    W_crew1=W_crew1*0.4536*g # N
    W_crew2=np.linspace(0,400, 20) # lbs
    W_crew2=W_crew2*0.4536*g # N 
    for i in range(len(W_crew1)):
        for j in range(len(W_crew2)):
            W_crew=W_crew1[i]+W_crew2[j]

            # On calcule le x_cg sans ballasts et on vérifie la stabilité
            W_b=0
            x_b=0
            #x_cg=compute_x_cg(W_b, x_b, W_crew1[i], W_crew2[j])
            #h=compute_h(x_cg)
            #h_n=compute_hn(x_cg)
            #K_n,is_kn_ok=compute_K_n(h_n,h)
            K_n, is_kn_ok = compute_K_n_bis(W_b, x_b, W_crew1[i], W_crew2[j])

            if is_kn_ok==True:
                plt.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color='#7db928')

            if is_kn_ok==False:
              m_lim_tot=60 #Total ballast weight admissible
              m_lim_tail=15 #Maximum ballast weight at the tail 
              m_lim_crew=m_lim_tot-m_lim_tail #Maximum ballast weight under the crew
              mb1=np.arange(0,m_lim_crew+1,7.5)*g
              mb2=np.arange(0,m_lim_tail+1,7.5)*g
              fine=False
              color="hello"
              for k in range(len(mb1)):
                  if(fine==True):
                      break
                  for l in range(len(mb2)):
                      W_b=mb1[k]+mb2[l]
                      if(W_b==0):
                          continue
                      x_b=(x_b_w*mb1[k]+x_b_t*mb2[l])/W_b # center of gravity of the ballasts 
                      #x_cg=compute_x_cg(W_b, x_b, W_crew1[i], W_crew2[j])
                      #h=compute_h(x_cg)
                      #h_n=compute_hn(x_cg)
                      #K_n,is_kn_ok=compute_K_n(h_n,h)
                      K_n, is_kn_ok = compute_K_n_bis(W_b, x_b, W_crew1[i], W_crew2[j])
                      if(is_kn_ok==True):
                          if(mb1[k]==0):
                              color='orange'
                          if(mb2[l]==0):
                              color='deepskyblue'
                          else:
                              color='orange'
                          fine=True
                          break
                
              if(fine==False):
                  
                  plt.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color='#e62d31')
              else:
                  plt.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color=color)
    red_patch = mpatches.Patch(color='#e62d31', label='Unstable')
    green_patch = mpatches.Patch(color='#7db928', label='Stable')
    orange_patch = mpatches.Patch(color='orange', label='Tail Ballast')
    blue_patch = mpatches.Patch(color='deepskyblue',label='Crew Ballast')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17),ncol=3, fancybox=True, shadow=True,fontsize=11,handles=[red_patch,green_patch,orange_patch,blue_patch])                            
    plt.xlabel("Crew member 2 [lb]",fontsize=11)
    plt.ylabel("Crew member 1 [lb]",fontsize=11)
    plt.savefig("cgdiag.pdf",bbox_inches='tight')
    plt.show()
    return 

# Stability diagram
cg_diag()



