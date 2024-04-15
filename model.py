import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator, interp2d
from scipy.optimize import fmin
from data import *
from airfoil import *
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

class Wing:
    @property
    def taper(self):
        return self.c_tip/self.c_root
    
    @property
    def S(self):
        return (self.c_tip + self.c_root) * self.b / 2
    
    @property
    def AR(self):
        return self.b**2 / self.S

    @property
    def ARh(self):
        return self.AR * (1 + 2*self.wl_h/self.b)**2

    @property
    def c_mac(self):
        return 2/3 * self.c_root*(1 + self.taper + self.taper**2)/(1 + self.taper)

    def Re_mac(self, ctx):
        return self.c_mac * ctx.V * ctx.rho / ctx.mu

    def sweep(self, x):
        x_tip = x*self.c_tip + self.tip_offset
        dx_root = x*self.c_root - x_tip
        return -np.arctan(dx_root / self.b)

    def CL_alpha(self, ctx):
        return self.mCL_alpha(self, ctx)

    def CL(self, alpha, ctx):
        return self.mCL(self, alpha, ctx)

    def CDv(self, alpha, ctx):
        return self.mCDv(self, alpha, ctx)
    
    def CDif(self, alpha, ctx):
        return self.mCDif(self, alpha, ctx)
    
    def CDin(self, alpha, ctx):
        return self.mCDin(self, alpha, ctx)

    def CD(self, alpha, ctx):
        return self.mCDv(self, alpha, ctx) + self.mCDif(self, alpha, ctx) + self.mCDin(self, alpha, ctx)

def WingCL_alpha_DATCOM(w, ctx):
    M = ctx.V / 330
    beta = np.sqrt(1 - M**2)
    k = w.foil.Cl_alpha(ctx.alpha_e, w.Re_mac(ctx)) / (2*np.pi/beta)
    return (2*np.pi*w.AR)/(2 + np.sqrt((w.AR*beta/k)**2 * (1 + np.tan(w.sweep(0.5))**2 / beta**2) + 4))

def WingCL_alpha_Anderson(w, ctx):
    a0 = w.foil.Cl_alpha(ctx.alpha_e, w.Re_mac(ctx))
    return 0.995 * a0 / (1 + 2*w.taper/w.AR/(1 + w.taper) + a0/np.pi/w.AR)

def WingCL_DATCOM(w, alpha, ctx):
    CL_alpha = w.CL_alpha(ctx)
    alpha0 = w.foil.alpha0(w.Re_mac(ctx))

    return CL_alpha * (alpha - alpha0)

def WingCDif_Zero(w, alpha, ctx):
    return 0

class WingCDin_Oswald:
    def __init__(self, e):
        self.e = e

    def __call__(self, w, alpha, ctx):
        CL = w.CL(alpha, ctx)
        return CL**2 / (np.pi * self.e * w.AR)

class WingCDv_Strip:
    def __init__(self, n, half=True):
        self.n = n
        self.half = half

    def __call__(self, w, alpha, ctx):
        d = 2 if self.half else 1

        # Calculate wing drag
        x = np.linspace(0, w.b/d, self.n)
        xr = np.linspace(0, 1, self.n)
        dx = x[1:] - x[:-1]
        cn_root = (w.c_root * (1 - xr[:-1]) + w.c_tip * xr[:-1])
        cn_tip = (w.c_root * (1 - xr[1:]) + w.c_tip * xr[1:])
        cn_mid = (cn_root + cn_tip)/2
        dS = (cn_root + cn_tip) * dx / 2
        Re = cn_mid * ctx.V * ctx.rho / ctx.mu
        Cds = w.foil.Cd(alpha, Re)*dS

        # Calculate wingtip drag
        S_wl = (w.c_tip + w.wl_tip) * w.wl_h / 2
        Cd_wl = 2 * w.foil.Cd(ctx.beta, Re[-1]) * S_wl / w.S

        return np.sum(Cds)/(w.S/d) + Cd_wl

fx62k131 = Airfoil(
    {
        0.100e6: 'FoilData/FX62K131/T1_Re0.100_M0.00_N5.0.csv',
        0.130e6: 'FoilData/FX62K131/T1_Re0.130_M0.00_N5.0.csv',
        0.160e6: 'FoilData/FX62K131/T1_Re0.160_M0.00_N5.0.csv',
        0.330e6: 'FoilData/FX62K131/T1_Re0.330_M0.00_N5.0.csv',
        0.500e6: 'FoilData/FX62K131/T1_Re0.500_M0.00_N5.0.csv',
        0.750e6: 'FoilData/FX62K131/T1_Re0.750_M0.00_N5.0.csv',
        1.000e6: 'FoilData/FX62K131/T1_Re1.000_M0.00_N5.0.csv',
        1.500e6: 'FoilData/FX62K131/T1_Re1.500_M0.00_N5.0.csv',
        1.750e6: 'FoilData/FX62K131/T1_Re1.750_M0.00_N5.0.csv',
        3.000e6: 'FoilData/FX62K131/T1_Re3.000_M0.00_N5.0.csv',
    }
)

naca0012 = Airfoil(
    {
        0.100e6: 'FoilData/NACA0012/T1_Re0.100_M0.00_N5.0.csv',
        0.130e6: 'FoilData/NACA0012/T1_Re0.130_M0.00_N5.0.csv',
        0.160e6: 'FoilData/NACA0012/T1_Re0.160_M0.00_N5.0.csv',
        0.330e6: 'FoilData/NACA0012/T1_Re0.330_M0.00_N5.0.csv',
        0.500e6: 'FoilData/NACA0012/T1_Re0.500_M0.00_N5.0.csv',
        0.750e6: 'FoilData/NACA0012/T1_Re0.750_M0.00_N5.0.csv',
        1.000e6: 'FoilData/NACA0012/T1_Re1.000_M0.00_N5.0.csv',
        1.500e6: 'FoilData/NACA0012/T1_Re1.500_M0.00_N5.0.csv',
        1.750e6: 'FoilData/NACA0012/T1_Re1.750_M0.00_N5.0.csv',
        3.000e6: 'FoilData/NACA0012/T1_Re3.000_M0.00_N5.0.csv',
    }
)


class Fuselage:
    def Re(self, ctx):
        return ctx.rho*ctx.V*self.L/ctx.mu

    def CDv(self, alpha, ctx):
        return self.mCDv(self, alpha, ctx)
    def CDif(self, alpha, ctx):
        return self.mCDif(self, alpha, ctx)
    def CDin(self, alpha, ctx):
        return self.mCDin(self, alpha, ctx)
    def CDmisc(self, alpha, ctx):
        return self.mCDmisc(self, alpha, ctx)

def Fuselage_CDin(f, alpha, ctx):
    CD_vortex = 0.15 * alpha**2 * f.Vf**(2/3) / w.S

    return CD_vortex

def Fuselage_CDv(f, alpha, ctx):
    dfeff = np.sqrt(4/np.pi * Ac)
    leff = l_fus / dfeff
    phi_f = 60/leff**3 + 0.0025*leff # DATCOM
    #phi_f = 2.2 / leff**1.5 + 3.8 / leff**3 # Torenbeek

    # Base drag
    Cf = 0.455 / np.log10(f.Re(ctx))**2.58
    CDb = Cf * S_fuselage * (1 + phi_f) / w.S

    # Due to angle of attack
    CDa = k2k1/w.S*2*alpha**2*S0 + 2*alpha**3/w.S/2/np.pi*d_eta*d_cdc*Sd

    return CDb + CDa

def Fuselage_CDif(f, alpha, ctx):
    # Wing wash interference
    nf = 0.637/w.b
    CL_o = w.CL(0, ctx)
    CD_wash = 0.55*nf/(1 + w.taper) * (2 - np.pi*nf)*CL_o**2 / np.pi / w.AR

    # Wing viscous interference
    Cf = 0.455 / np.log10(w.Re_mac(ctx))**2.58
    CD_viscous = 1.5 * Cf * 0.131 * 4.5 * 4.5 * w.c_root * np.cos(w.sweep(0.5))**2 / w.S
    
    # Tail wash interference
    CD_wash_tail = 0
    
    # Tail viscous interference
    CD_viscous_tail = 0.1*(h.CDv(alpha, ctx)*h.S/w.S + v.CDv(alpha, ctx)*v.S/w.S)
    
    return CD_wash + CD_viscous + CD_wash_tail + CD_viscous_tail

def Fuselage_CDmisc(f, alpha, ctx):
    # Wheels
    fr = 0.3556*0.11938 + 0.1524*0.03048
    CD_wheels = 0.22*fr/w.S

    # Cockpit/canopy (well streamlined cylindrical mid-section)
    S_canopy = 0.4*Ac # approximation
    phi_f = 0
    CD_canopy = 0.85*(1 + phi_f)*(1 + 4.5*alpha)*(0.075 + 0.02)*S_canopy/w.S
    
    # Imperfections
    CD_wing = 0.06*w.CDv(alpha, ctx)
    CD_fus = 0.07*fus.CDv(alpha, ctx)

    # Retracted controls
    CD_ctrl_wing = 1.8e-4
    CD_ctrl_fin = 1.8e-4*v.S/w.S
    CD_ctrl_elevator = 1.8e-4*h.S/w.S
    return CD_wheels + CD_canopy + CD_wing + CD_fus + CD_ctrl_wing + CD_ctrl_fin + CD_ctrl_elevator

# t = np.linspace(0, 10, 100)
# v = 1 + 2.3*t**2
# plt.plot(t, v)
# plt.fill_between(t, v, alpha=0.5)
# pv = v.copy()
# v += 1.04*t**2
# plt.plot(t, v)
# plt.fill_between(t, v, pv, alpha=0.5)
# pv = v.copy()
# v += 20
# plt.plot(t, v)
# plt.fill_between(t, v, pv, alpha=0.5)
# plt.show()

class Ctx:
    pass

ctx = Ctx()
ctx.alpha = 0/180 * np.pi
ctx.alpha_e = 0
ctx.beta = 0
ctx.V = 32
ctx.rho = 1.2
ctx.mu = 1.7e-5

v = Wing()
v.foil = naca0012
v.c_tip = c_v_tip
v.c_root = c_v_root
v.b = b_v
v.wl_h = 0
v.wl_tip = 0
v.mCL_alpha = WingCL_alpha_DATCOM
v.mCL = WingCDif_Zero
v.mCDv = WingCDv_Strip(50, False)
v.mCDif = WingCDif_Zero
v.mCDin = WingCDif_Zero

h = Wing()
h.foil = naca0012
h.c_tip = c_h_tip
h.c_root = c_h_root
h.b = b_h
h.wl_h = 0
h.wl_tip = 0
h.tip_offset = 0
h.mCL_alpha = WingCL_alpha_DATCOM
h.mCL = WingCDif_Zero
h.mCDv = WingCDv_Strip(50)
h.mCDif = WingCDif_Zero
h.mCDin = WingCDin_Oswald(0.9)

w = Wing()
w.foil = fx62k131
w.c_tip = c_w_tip
w.c_root = c_w_root
w.b = b_w
w.wl_h = 0.5
w.wl_tip = w.c_tip * w.taper
w.tip_offset = 0
w.mCL_alpha = WingCL_alpha_Anderson
w.mCL = WingCL_DATCOM
w.mCDv = WingCDv_Strip(50)
w.mCDif = WingCDif_Zero
w.mCDin = WingCDin_Oswald(0.9335)


fus = Fuselage()
fus.Vf = Vf
fus.L = l_fus
fus.mCDin = Fuselage_CDin
fus.mCDv = Fuselage_CDv
fus.mCDif = Fuselage_CDif
fus.mCDmisc = Fuselage_CDmisc

@np.vectorize
def CDv(alpha, ctx):
    a = alpha
    aw = a + theta_w
    ah = a - theta_h
    return w.CDv(aw, ctx) + h.CDv(ah, ctx)*h.S/w.S + v.CDv(0, ctx)*v.S/w.S + fus.CDv(a, ctx)
@np.vectorize
def CDin(alpha, ctx):
    a = alpha
    aw = a + theta_w
    ah = a - theta_h
    return w.CDin(aw, ctx) + h.CDin(ah, ctx)*h.S/w.S + v.CDin(0, ctx)*v.S/w.S + fus.CDin(a, ctx)
@np.vectorize
def CDif(alpha, ctx):
    a = alpha
    aw = a + theta_w
    ah = a - theta_h
    return w.CDif(aw, ctx) + h.CDif(ah, ctx)*h.S/w.S + v.CDif(0, ctx)*v.S/w.S + fus.CDif(a, ctx)
@np.vectorize
def CDmisc(alpha, ctx):
    a = alpha
    aw = a + theta_w
    ah = a - theta_h
    return fus.CDmisc(a, ctx)

def CD(alpha, ctx):

    return w.CD(alpha, ctx) + h.CD(alpha - 1.8/180*np.pi, ctx)*h.S/w.S + v.CD(0, ctx)*v.S/w.S + fus.CDv(alpha, ctx) + fus.CDin(alpha, ctx) + fus.CDif(alpha, ctx) + fus.CDmisc(alpha, ctx)

def CD2(alpha, ctx):
    return CDv(alpha, ctx) + CDin(alpha, ctx) + CDif(alpha, ctx) + CDmisc(alpha, ctx)

a = 0/180*np.pi
ah = a - 1.8/180*np.pi
print('CD', CD(0/180*np.pi, ctx))
print('CD2', CD2(0/180*np.pi, ctx))

print('RE', fus.Re(ctx))


alpha = np.linspace(-1, 5, 100)/180*np.pi
cdv = CDv(alpha, ctx)
plt.plot(alpha/np.pi*180, cdv, label='$C_{D,\\mathrm{v}}$')
plt.fill_between(alpha/np.pi*180, cdv, alpha=0.5)
cdin = cdv + CDin(alpha, ctx)
plt.plot(alpha/np.pi*180, cdin, label='$+ C_{D,\\mathrm{in}}$')
plt.fill_between(alpha/np.pi*180, cdin, cdv, alpha=0.5)
cdif = cdin + CDif(alpha, ctx)
plt.plot(alpha/np.pi*180, cdif, label='$+ C_{D,\\mathrm{if}}$')
plt.fill_between(alpha/np.pi*180, cdif, cdin, alpha=0.5)
cdmisc = cdif + CDmisc(alpha, ctx)
plt.plot(alpha/np.pi*180, cdmisc, label='$+ C_{D,\\mathrm{misc}}$')
plt.fill_between(alpha/np.pi*180, cdmisc, cdif, alpha=0.5)
plt.xlabel('Angle of attack [°]')
plt.ylabel('Drag coefficient $C_D$')
plt.legend()
plt.show()


exit(0)


# print('Re_stall', w.Re_mac(ctx))
# print('Cl_max_stall', w.foil.Cl_max(w.Re_mac(ctx)))

# def downwash():
#     K_A = 1/AR_w - 1/(1+AR_w**1.7)
#     K_lambda = (10-3*lambda_w)/7
#     K_H = (1-h_H/b_w)/(np.sqrt(2*l_H/b_w))
#     grad_down_M0 = 4.44*(K_A*K_lambda*K_H*np.cos(w.sweep(0.25)))**1.19
#     grad_down = grad_down_M0*1      # approx that CL_alpha,wM/CL_alpha,wM=0 ~ 1
#     return grad_down 

# print('downwash', downwash())

# def get_KN():
#     CLw_alpha = w.CL_alpha(ctx)

#     return 2 * S_N / (CLw_alpha * w.S)

# SW_body = pe.b_i((x_debut_wing + c_w_root/2)*1000)/1000 * c_w_root
# print('SW_body', SW_body)

# w_coef = get_KN() + 1.03 + 0.05
# print('KN', get_KN())
# print('d/b wing', 0.8221144372/w.b)
# print('d/b tail', 0.8221144372/h.b)

# print((w.S - SW_body)/w.S)

# def CL_alpha():
#     w_alpha0 = w.foil.alpha0(w.Re_mac(ctx))
#     h_alpha0 = h.foil.alpha0(h.Re_mac(ctx))

#     return (w.S - SW_body)/w.S * w_coef * w.CL_alpha(ctx) + (1 - downwash()) * h.CL_alpha(ctx) * h.S / w.S


# def CL(alpha):
#     w_alpha0 = w.foil.alpha0(w.Re_mac(ctx))
#     h_alpha0 = h.foil.alpha0(h.Re_mac(ctx))

#     return (w.S - SW_body)/w.S * w_coef * w.CL_alpha(ctx) * (alpha + theta_w - w_alpha0) + (1 - downwash()) * h.CL_alpha(ctx) * h.S / w.S * (alpha + theta_h)

# def alpha0():
#     return fmin(lambda alpha: CL(alpha)**2, -5)[0]

# # V = np.linspace(10, 70, 100)
# # Re = np.zeros_like(V)
# # CL_alphas = np.zeros_like(V)
# # for i, v in enumerate(V):
# #     ctx.V = v
# #     Re[i] = w.Re_mac(ctx)
# #     CL_alphas[i] = CL_alpha()

# pdata = np.genfromtxt('T1-32_0 m_s-VLM2.csv', delimiter=',')

# alpha = np.linspace(-8, 8, 10) / 180 * np.pi
# plt.figure(figsize=(5, 2.5))
# plt.plot(alpha / np.pi * 180, CL(alpha), label='Empirical w/ fus.')
# w_coef = 1 #get_KN() + 1.03 + 0.05
# plt.plot(alpha / np.pi * 180, CL(alpha), label='Empirical w/o fus.')
# plt.plot(pdata[:,0], pdata[:,2], label='VLM')
# plt.xlabel('Angle of attack [°]')
# plt.ylabel('$C_L$')
# plt.text(-7, 0.732-0.669 + 0.575, 'Cruise $C_L$')
# plt.axhline(0.575, linestyle='--', c='gray')
# plt.legend()
# plt.tight_layout()
# plt.savefig('lift_cruise.pdf')
# plt.show()


# # Find the zero lift angle
# a = (pdata[5,2] - pdata[3,2])/((pdata[5,0] - pdata[3,0]))
# b = pdata[5,0] - pdata[5,2] / a
# print(a, b)

# print('empirical CL_alpha', CL_alpha(), alpha0() * 180 / np.pi)
# print('VLM CL_alpha', a / np.pi * 180, b)

# # plt.plot(Re, CL_alphas)
# # plt.show()


# #print(w.CD(ctx) + v.CD(ctx)*v.S/w.S + h.CD(ctx)*h.S/w.S)
# def polar():
#     alpha = np.linspace(-8, 8, 30)/180 * np.pi
#     CL = np.zeros_like(alpha)
#     CD = np.zeros_like(alpha)
#     for i, a in enumerate(alpha):
#         ctx.alpha = a
#         CL[i] = w.CL(ctx)
#         CD[i] = w.CD(ctx)

#     plt.plot(CD, CL)
#     plt.show()