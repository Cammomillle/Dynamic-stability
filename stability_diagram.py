import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.patches as mpatches
from data import *
from functions import *
from cycler import cycler

plt.rcdefaults()  # restore the defaults

custom_colors = ['#00707f', '#f07f3c', '#7db928', '#ffd000', '#e62d31', '#5b57a2']  # turquoise, orange, green, yellow, red, blue 
custom_cycler = cycler(color=custom_colors)
plt.rcParams.update({
    'axes.prop_cycle': custom_cycler,
    'text.usetex': True,
    'text.latex.preamble': '\\usepackage{stix2}',
    'font.family': 'serif',
    'font.serif': 'stix2'
})
#************************************************
#*** Pitch stability diagram of the sailplane *** 
#************************************************

def plot_polygone(ax,x_ok,y_ok,color):
    y_ok_unique=np.unique(y_ok)
    x_surf=[]
    y_surf=[]
    y_act=0
    x_temp_tab=[]
    x_elems=[]
    for i in range(len(x_ok)):
        if(i==0):
            continue
        if x_ok[i]>x_ok[i-1]:
            x_temp_tab.append(x_ok[i-1])
        else:
            x_temp_tab.append(x_ok[i-1])
            x_elems.append(np.copy(x_temp_tab))
            x_temp_tab=[]
    x_temp_tab.append(x_ok[-1])
    x_elems.append(x_temp_tab)
    x_surface=[]
    y_surface=[]
    x_surface_left=[]
    y_surface_left=[]
    y=np.unique(y_ok)
    for i in range(len(x_elems)):
          x_el=x_elems[i]
          x_surface.append(x_el[-1])
          y_surface.append(y[i])
          x_surface_left.append(x_el[0])
          y_surface_left.append(y[i])
    x_surface_left.reverse()
    y_surface_left.reverse()
    x_surface=x_surface+x_surface_left
    y_surface=y_surface+y_surface_left
    x_surface.append(x_surface[0])
    y_surface.append(y_surface[0])
    ax.fill(x_surface,y_surface,facecolor=color,edgecolor='black',linewidth=0.3,alpha=1)
def cg_diag():
    W_crew1=np.linspace(0,400,20) # lbs
    W_crew1=W_crew1*0.4536*g # N
    W_crew2=np.linspace(0,400,20) # lbs
    W_crew2=W_crew2*0.4536*g # N 
    x_ok=[]
    x_ok_ballast_tail=[]
    x_ok_ballast_crew=[]
    y_ok=[]
    y_ok_ballast_tail=[]
    y_ok_ballast_crew=[]
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    for i in range(len(W_crew1)):
        print("% done", W_crew1[i]/W_crew1[-1])
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
                ax2.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color='#7db928')
                x_ok.append(W_crew2[j]/(0.4536*g))
                y_ok.append(W_crew1[i]/(0.4536*g))

            if is_kn_ok==False:
              m_lim_tot=45 #Total ballast weight admissible
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
                              color='orange'
                          else:
                              color='orange'
                          fine=True
                          break
                
              if(fine==False):
                  ax2.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color='#e62d31')
              else:
                  ax2.scatter(W_crew2[j]/(0.4536*g),W_crew1[i]/(0.4536*g),color=color)
                  x_ok.append(W_crew2[j]/(0.4536*g))
                  y_ok.append(W_crew1[i]/(0.4536*g))
                  """if(color=='deepskyblue'):
                    x_ok_ballast_tail.append(W_crew2[j]/(0.4536*g))
                    y_ok_ballast_tail.append(W_crew1[i]/(0.4536*g))
                  else:
                    x_ok_ballast_crew.append(W_crew2[j]/(0.4536*g))
                    y_ok_ballast_crew.append(W_crew1[i]/(0.4536*g))"""
    x=[0,400,400,0,0]
    y=[0,0,400,400,0]
    ax1.fill(x,y,facecolor='#e62d31',edgecolor='black',linewidth=0.3,alpha=1)
    plot_polygone(ax1,x_ok,y_ok,color='#7db928')
    ax1.grid(visible=True,color='black', linestyle='-', linewidth=0.5)
    #plot_polygone(ax1,x_ok_ballast_crew,y_ok_ballast_crew,color='orange')
    #plot_polygone(ax1,x_ok_ballast_tail,y_ok_ballast_tail,color='orange')
    red_patch = mpatches.Patch(color='#e62d31', label='Unstable')
    green_patch = mpatches.Patch(color='#7db928', label='Stable')
    orange_patch = mpatches.Patch(color='orange', label='Stable with tail ballast')
    blue_patch = mpatches.Patch(color='deepskyblue',label='Stable with crew ballast')
    ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17),ncol=3, fancybox=True, shadow=True,fontsize=11,handles=[red_patch,green_patch,orange_patch])                            
    ax2.set_xlabel("Crew member 2 [lb]",fontsize=11)
    ax2.set_ylabel("Crew member 1 [lb]",fontsize=11)
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.17),ncol=2, fancybox=True, shadow=True,fontsize=11,handles=[red_patch,green_patch])                            
    ax1.set_xlabel("Crew member 2 [lb]",fontsize=11)
    ax1.set_ylabel("Crew member 1 [lb]",fontsize=11)
    fig1.savefig("cgdiag.pdf",bbox_inches='tight')
    plt.show()
    return 

# Stability diagram
cg_diag()



