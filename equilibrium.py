from data import *
from functions import *
#from functions_old import * #a décommenter pour voir avec l'ancienne version
from matplotlib import pyplot as plt
#from stability_test import ballast_positions

plt.rcdefaults()  # restore the defaults
from cycler import cycler


custom_colors = ['#00707f', '#f07f3c', '#7db928', '#ffd000', '#e62d31', '#5b57a2']
custom_cycler = cycler(color=custom_colors)
plt.rcParams.update({
    'axes.prop_cycle': custom_cycler,
    'text.usetex': True,
    'text.latex.preamble': '\\usepackage{stix2} \\usepackage{GoSans}',
    'font.family': 'sans-serif',
    'font.serif': 'GoSans'
})

#***********************************************************************************
#*** Equilibrium for a configuration (number of pilot, pilots' weight, ballasts) *** 
#***********************************************************************************

#******* Weights [kilograms * g] **********
W_crew1 = 80*g   # 1 st crew
W_crew2 = 80*g   # 2 nd crew

#compute the ballasts' positions and weights given the weight of the passenger
def ballast_positions(x_crew1,x_crew2, x_b_w, x_b_t, W_crew1, W_crew2, print_T=True):
    w_b_1=w_b_2=0
    W_crew=W_crew1+W_crew2
    x_crew=0
    if(W_crew!=0):
        x_crew=(W_crew1*x_crew1+W_crew2*x_crew2)/W_crew

    # On calcule le x_cg sans ballasts et on vérifie la stabilité
    W_b=0
    x_b=0
    x_cg=compute_x_cg(W_b, x_b, W_crew1, W_crew2)
    h=compute_h(x_cg)
    h_n=compute_hn(x_cg)
    K_n,is_kn_ok=compute_K_n(h_n,h)
    #print(x_cg, h, h_n)

    if is_kn_ok==True and print_T:
        print("Stable without ballasts")

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
              x_cg=compute_x_cg(W_b, x_b, W_crew1, W_crew2)
              h=compute_h(x_cg)
              h_n=compute_hn(x_cg)
              K_n,is_kn_ok=compute_K_n(h_n,h)
              if(is_kn_ok==True):
                  fine=True
                  W_b_list.append(W_b)
                  best_list.append([mb1[k],mb2[l]])

        try:          
            Wbmin=W_b_list[np.argmin(W_b_list)]
            print("[weight b1, weight b2, CG with ballasts]: ", best_list)
        except:
            print("No stable configurations with ballasts")
            return x_cg
        for val in best_list:
          if print_T:print('ballast that we keep: ',val)
          if val[0]+val[1]==Wbmin:
              w_b_1=val[0]
              w_b_2=val[1]
              break
                  
        if(fine==False and print_T):
          print("No stable configurations with ballasts")
        elif print_T:
          print("Ballasts : {0} kg at the wing CG and ".format(w_b_1/g)+"{0} kg at 3/4 of the tail length".format(w_b_2/g))
    return w_b_1, w_b_2 

#compute the wing lift curve slope for a swept wing config (slide 121 L1 aéro)
def Wing_lift_curve_slope(AR, cl_alpha, Sweepback_angle_half, beta=1):
    k=beta*cl_alpha/(2*np.pi)
    Lambda_beta=Sweepback_angle_half
    return 2*np.pi/(beta*(2/(beta*AR)+np.sqrt(1/(k**2*(np.cos(Lambda_beta))**2)+(2/(beta*AR))**2)))
print("ssssss",Wing_lift_curve_slope(AR_h, dclh_alpha_h, sweep_h_half))
print("ssssss",Wing_lift_curve_slope(AR_v, dclF_dbeta, sweep_v))


#compute the lift of the horizontal stabilizer
#Different formulas. Because now the empennage is swept, we use C_L2 and L2
def Lift_horizontal_stabilizer(alpha):
    alpha=np.pi/180*alpha
    S_emp=S_h
    AR=AR_h
    
    #Wing maximum lift (empennage) coefficient (slide 122 L1 aéro )
    k_s=0.95            #tapered wing 
    c_l_max_root=1.3    #to verify !!!!! => ça vient d'où ?? wing ??? demander
    c_l_alpha=dclh_alpha_h
    
    C_L_max=k_s*c_l_max_root    
    L_max=0.5*rho*V0**2*S_emp*C_L_max
    
    #Wing lift (empennage) for unswept wing (slide 120 L1 aéro)
    E=1+2*lambda_h/(AR*(1+lambda_h))
    C_L_alpha=.995*c_l_alpha/(E+c_l_alpha/(np.pi*AR))   
    
    C_L=alpha*C_L_alpha
    L=0.5*rho*V0**2*S_emp*C_L
    
    #Wing lift (empennage) for swept wing (slide 121 L1 aéro)
    C_L2=alpha*Wing_lift_curve_slope(AR, c_l_alpha, sweep_h_half)
    L2=0.5*rho*V0**2*S_emp*C_L2
    
    #print("\n for alpha= ", alpha*180/np.pi,"C_L= ", C_L, "C_L_2= ", C_L2, "L_max= ", L_max)
    return L2, C_L2

#compute the drag of the horizontal stabilizer
def Drag_horizontal_stabilizer(alpha, C_L):
    alpha=np.pi/180*alpha
    e=0.971 # correction, (data d'un glider) dans le slide 61 du cours concptual design d'autres valeurs sont données
    AR=AR_h
    S_emp=S_h
    C_D_0=CD_0_h     # we take the same value as for the wing
    C_D=C_D_0+C_L**2/(e*np.pi*AR)   
    D_horiz=0.5*rho*V0**2*S_emp*C_D
    return D_horiz

#first version, not used for now
def compute_moments(alpha):
    #Compute the moment at the MAC for the wing
    N_wing=np.cos(aoa_w*np.pi/180)*L_w + np.sin(aoa_w*np.pi/180)*D_w #L=Lift, D=Drag
    #M_AC_wing=(x_AC-x_CP)*N_wing
    M_AC_wing=0     # ok comme assumption ?? on est en 3D ici
    
    #Compute the moment at the MAC for the wing
    alpha_horiz=alpha
    L_horiz, C_L_horiz=Lift_horizontal_stabilizer(alpha_horiz)
    D_horiz=Drag_horizontal_stabilizer(alpha_horiz, C_L_horiz)
    N_horiz=np.cos(alpha_horiz*np.pi/180)*L_horiz + np.sin(alpha_horiz*np.pi/180)*D_horiz #L=Lift, D=Drag
    #M_AC_horiz=(x_AC-x_CP)*N_horiz
    M_AC_horiz=0    # ok comme assumption ?? on est en 3D ici
    return M_AC_wing, M_AC_horiz

#formule d'alejandro ???
def compute_c_ac(taper, c_root=c_w_root):  
    c_ac = 2/3 * c_root * (1 + taper + taper**2)/(1 + taper)
    return c_ac

#compute moment
def compute_moments_2():
    #pour la wing, c_m ~ -0.18 (quarter chord)
    c_ac=compute_c_ac(lambda_w, c_w_root)
    M_wing_AC = 0.5*rho*V0**2*S_w_total*c_ac*c_m_w
    #print(M_wing_AC)
    
    #pour l'empennage, c_m ~ 0 (1ere approx, cf airfoil tools, Re=1.000.000) !!!(quarter chord)
    c_m=0 
    c_ac=compute_c_ac(lambda_h, c_h_root)
    #print("xxxx",c_ac)
    M_h_ac=0.5*rho*V0**2*S_h*c_ac*c_m
    #print(M_wing_AC, M_h_ac)
    return M_wing_AC, M_h_ac


#En fonction de l'angle alpha (degré) de l'horizontal stabilizer, renvoie le déséquilibre de force et de moment
def Equilibrium_validation(alpha, W_crew11, W_crew22, print_T=True):
    #Equilibre des forces
    #Sans ballast
    W_b = 0
    x_b = 0
    W_total=compute_total_mass(W_b, W_crew11, W_crew22)
        
    L_T,C=Lift_horizontal_stabilizer(alpha)
    #print(L_T)
    Delta_Force_wo_b=L_w+L_T-W_total
    
    x_cg=compute_x_cg(W_b, x_b, W_crew11, W_crew22)
    h=compute_h(x_cg)
    x_le_wing = compute_x_mac(x_debut_wing, b_w, lambda_w, sweep_w_le)
    h_0=(x_ac_w-x_le_wing)/c_mac_w
    
    l_T=x_ac_h-x_cg     
    M_0, M_T=compute_moments_2() #calcul les moments de la wing et de l'empennage (sens positif pour le nez de l'avion qui monte)
    print("M_0:", M_0)
    #print("h, h_0, c__w:", h, h_0, c__w)
    Delta_Moment_wo_b = M_0 + L_w*(h-h_0)*c_mac_w - L_T*l_T + M_T

    #Avec ballast
    W_ballast1, W_ballast2=ballast_positions(x_crew1,x_crew2, x_b_w, x_b_t, W_crew11, W_crew22, print_T)
    #W_ballast1=W_ballast1
    W_total=compute_total_mass((W_ballast1 + W_ballast2), W_crew11, W_crew22)
    print("W_tot: ", W_total)
    L_T,C=Lift_horizontal_stabilizer(alpha)
    print("L_T:", L_T, alpha)
    Delta_Force_w_b=L_w+L_T-W_total 
    L_w_bis=W_total-L_T # Autre approche: tiens compte de la variation de lift produite par la wing en fonction du chargement
    print("xxxx", L_w_bis)
    x_b = 0
    if(W_ballast1 + W_ballast2 != 0):
        x_b = (W_ballast1*x_b_w + W_ballast2*x_b_t)/(W_ballast1+W_ballast2)

    x_cg=compute_x_cg(W_ballast1+W_ballast2, x_b, W_crew11, W_crew22)
    h=compute_h(x_cg)
    l_T=x_ac_h-x_cg

    M_0, M_T=compute_moments_2() #calcul les moments de la wing et de l'empennage (sens positif pour le nez de l'avion qui monte)
    
    #print("h, h_0, c__w:", h, h_0, c__w)
    Delta_Moment_w_b=M_0+L_w_bis*(h-h_0)*c_mac_w-L_T*l_T+M_T
    if(print_T==True):
        print("Delta_Force without ballast : ", Delta_Force_wo_b)
        print("Delta_Force with ballast    : ", Delta_Force_w_b)
        print("Delta_Moment without ballast: ", Delta_Moment_wo_b)
        print("Delta_Moment with ballast   : ", Delta_Moment_w_b)
        print("x_CG :", x_cg)
    return Delta_Moment_w_b

#Equilibrium study
#Equilibrium_validation(0, 120*g,120*g)

#print(Lift_horizontal_stabiliser(5))
min_w=50
max_w=120
step=10

def Equilibrium_testing(alpha, print_T=False):
    tab=[]
    for i in range (min_w,max_w+1, step):
        if(print_T==True):print("============================\n 1 passenger weight=",i)
        tab.append(Equilibrium_validation(alpha, i*g, 0, print_T))
    for i in range (min_w,max_w+1, step):
        if(print_T==True):print("============================\n 2 passenger weight=",i)
        tab.append(Equilibrium_validation(alpha, i*g, i*g, print_T))
    lift_tail=Lift_horizontal_stabilizer(alpha)
    return tab, lift_tail
print_T=False
print(print_T)
alpha_tab=[0, -2, -3, -4]
t1, lift_tail_1=Equilibrium_testing(alpha_tab[0], print_T)
t2, lift_tail_2=Equilibrium_testing(alpha_tab[1], print_T)
t3, lift_tail_3=Equilibrium_testing(alpha_tab[2], print_T)
t4, lift_tail_4=Equilibrium_testing(alpha_tab[3], print_T)
size= int(len(t1)/2)
t00=np.linspace(min_w, max_w,size)
t01=np.zeros(size)
t02=np.hstack((t01, t00))
t03=np.hstack((t00, t00))

t=np.vstack((t03, t02, np.array(t1), np.array(t2), np.array(t3), np.array(t4)))
t=np.transpose(t)


#Analyse des résultats
def plot_moment_imbalance(t):
    for i in range(2, np.size(t[0]), 1):
        mid=int(len(t[:,i])/2)
        
        abscisse=np.arange(0, mid, 1)
        plt.plot((t[abscisse,0])*2.20462, (t[abscisse,i])*0.7376, '-o', label=r'$\rm \theta_h=$'+str(alpha_tab[i-2])+'°, 1 pass.', linewidth=0.5)
        plt.plot((t[abscisse+mid,0])*2.20462, (t[abscisse+mid,i])*0.7376, '-x', label=r'$\rm \theta_h=$'+str(alpha_tab[i-2])+'°, 2 pass.', color=plt.gca().lines[-1].get_color(), linewidth=0.5)
        plt.legend(ncol=2, columnspacing=0.3, loc='best')
        plt.xlabel("Weight of the crew members [lb]")
        plt.ylabel("Moment imbalance [ft-lb]")
plot_moment_imbalance(t)
plt.savefig("Moment_imbalance.pdf", bbox_inches='tight')    
plt.show()               
def print_table_to_latex(t):
    for i in range(len(t)):
        for j in range (np.size(t[0])):
            if (j!=(np.size(t[0])-1)):
                print(int(t[i, j]), "&", end="")
            else:
                print(int(t[i, j]), end="")
        print("\\\\")
#print_table_to_latex(t)
print(lift_tail_1, lift_tail_2, lift_tail_3, lift_tail_4)


print("zzzzzzzzzzzzzzzzzzzzzzzzzz")
Equilibrium_validation(0, 80*g, 80*g)

print(Lift_horizontal_stabilizer(-2))


