import numpy as np

hc = 8.5
x1 = 3.5015
x0 = hc*(0.378 + 0.527*(x1/hc))

Amaxa = 1.08373/2
Amaxb = 0.69862/2
Amax = np.pi*Amaxa*Amaxb

S0a = 0.35311/2
S0b = 0.1005
S0 = np.pi*S0a*S0b

k2k1 = 0.8296470588235295
Vb = 1575323985.07004523 / 1e9

d_eta = 0.6826923076923079
d_cdc = 1.2

C0 = (0.89 + 0.64)/2
Sd = (hc-x0)*C0
S_fus = 12

Sw = 18
Swl = 0.9
ARe = 24.50000

f = hc/np.sqrt(4/np.pi * Amax)
IF_fus = 1.1
FF_fus = 1 + 60/f**3 + 0.0025*f

rho = 1.11165
V = 32
mu = 1.7575071729834082e-05
Re_fus = rho*V*hc/mu

Cf_fus = 0.455 / np.log10(Re_fus)**2.58

Sh = 1.80000
St = 1.58000
ARht = 8.88889
e = 0.895

def foildrag(aoa):
    x = aoa/np.pi*180
    a = 2.10572888e-09
    b = 5.46557733e-08
    c = 1.97321338e-04
    d = 5.58636985e-03
    return a*x**6 + b*x**4 + c*x**2 + d

def htfoildrag(aoa):
    x = aoa/np.pi*180
    a = 4.41122032e-09
    b = -1.23974674e-06 
    c = 1.74486543e-04
    d = 4.93028742e-03
    return a*x**6 + b*x**4 + c*x**2 + d

def wing_aoa(aoa):
    theta_wing = 1.17061/180 * np.pi
    return aoa + theta_wing

def htail_aoa(aoa):
    theta_wing = -2.64053/180 * np.pi
    return aoa + theta_wing

# Fuselage Cl (ref = Sw; assumption)
def Cl_fus(aoa):
    return 0

# Fuselage Cd (ref = Sw; DATCOM)
def Cd_fus(aoa):
    return Cf_fus*IF_fus*FF_fus*S_fus/Sw + k2k1/Sw*2*aoa**2*S0 + 2*aoa**3/Sw/2/np.pi*d_eta*d_cdc*Sd

# Wing Cl (ref = Sw; Anderson)
def Cl_wing(aoa):
    alpha0 = -5.88608/180 * np.pi
    a = 5.63912
    return a*(aoa - alpha0)

# Wing Cd (ref = Sw; X-FOIL)
def Cd_wing(aoa):
    IFw = 1.02
    IFwl = 1.01
    return IFw*(Cl_wing(aoa)**2 / (e*np.pi*ARe) + foildrag(aoa)) + IFwl*foildrag(0)*Swl/Sw

# Wing Cm (ref = Sw c_mac; assumption)
def Cm_wing(aoa):
    return -0.18

# Horizontal tail Cl (ref = Sw; Anderson)
def Cl_htail(aoa):
    a = 5.06900
    return a*htail_aoa(aoa)*Sh/Sw

# Horizontal tail Cd (ref = Sw; X-FOIL)
def Cd_htail(aoa):
    IFht = 1.04
    r = IFht * (Cl_htail(aoa)**2 / (e * np.pi * ARht) + htfoildrag(aoa))
    return Sh/Sw * r

# Horizontal tail Cm (ref = Sw c_mac; assumption)
def Cm_htail(aoa):
    return 0

# Vertical tail Cd (ref = Sw; X-FOIL)
def Cd_vtail(aoa):
    IFvt = 1.04
    r = IFvt * htfoildrag(0)
    return St/Sw * r

# Structure Cl (ref = Sw)
def Cl(aoa):
    return Cl_wing(wing_aoa(aoa)) + Cl_htail(htail_aoa(aoa))

# Structure Cd (ref = Sw; Torenbeek/DATCOM)
def Cd(aoa):
    return Cd_wing(wing_aoa(aoa)) + Cd_htail(htail_aoa(aoa)) + Cd_vtail(aoa) + Cd_fus(aoa) + IF_fus*0.18*0.03*0.075/Sw

# Structure Cm (ref = Sw c_mac)
def Cm(aoa):
    # These lengths are normalized by total length hc
    h_wing = 0.34483
    h_tail = 1.00000
    h_center = 0.35287
    return Cm_wing(aoa) + Cm_htail(aoa) - (h_wing - h_center)*hc*Cl_wing(aoa) - (h_tail - h_center)*hc*Cl_htail(aoa)
