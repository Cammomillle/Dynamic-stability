import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator, interp2d
from scipy.optimize import fmin
from data import *
from airfoil import *
from cycler import cycler
import atmosphere as atm

custom_colors = ['#00707f', '#f07f3c', '#7db928', '#ffd000', '#e62d31', '#5b57a2']
custom_cycler = cycler(color=custom_colors)
plt.rcParams.update({
    'axes.prop_cycle': custom_cycler,
    'text.usetex': True,
    'text.latex.preamble': '\\usepackage{stix2,siunitx}',
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
    #print('Back', 0.22*0.1524*0.03048/w.S)
    #print('Rear', 0.22*0.3556*0.11938/w.S)

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
ctx.H = 1066.8
ctx.rho = atm.density(ctx.H)
ctx.mu = atm.viscosity(ctx.H)

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
w.b = b_w*0.9
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
def CL_VLM(alpha, ctx):
    a = 6.02
    alpha0 = -5.46/180*np.pi
    return a * (alpha - alpha0)

def CL_max(ctx):
    Re_root = ctx.V*ctx.rho*w.c_root/ctx.mu
    Re_tip = ctx.V*ctx.rho*w.c_tip/ctx.mu
    return 0.95*(w.foil.Cl_max(Re_tip) + w.foil.Cl_max(Re_root))/2

def CL_stall(m, ctx):
    V = ctx.V
    for _ in range(10):
        Re_root = V*ctx.rho*w.c_root/ctx.mu
        Re_tip = V*ctx.rho*w.c_tip/ctx.mu
        Cl_tip = w.foil.Cl_max(Re_tip)
        Cl_root = w.foil.Cl_max(Re_root)
        Cl_max = 0.95 * (Cl_tip + Cl_root)/2
        V = np.sqrt(m*9.81/(0.5*ctx.rho*w.S*Cl_max))
    return Cl_max, V

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
def CDin_test(alpha, ctx):
    a = alpha
    aw = a + theta_w
    ah = a - theta_h
    return w.CDin(aw, ctx) + h.CDin(ah, ctx)*h.S/w.S + v.CDin(0, ctx)*v.S/w.S
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
@np.vectorize
def CD0(alpha, ctx):
    return CDv(alpha, ctx) + CDif(alpha, ctx) + CDmisc(alpha, ctx)
@np.vectorize
def CD(alpha, ctx):
    return CDv(alpha, ctx) + CDin(alpha, ctx) + CDif(alpha, ctx) + CDmisc(alpha, ctx)


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

# w_coef = 1 #get_KN() + 1.03 + 0.05
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

#     return (w.S - SW_body)/w.S * w_coef * w.CL_alpha(ctx) * (alpha + theta_w - w_alpha0), (1 - downwash()) * h.CL_alpha(ctx) * h.S / w.S * (alpha + theta_h)

print('alpha0', w.foil.alpha0(w.Re_mac(ctx))/np.pi*180)
print('alpha_max', w.foil.alpha_max(w.Re_mac(ctx))/np.pi*180)
print('CL_max', CL_max(ctx))
print('a', w.CL_alpha(ctx))
print('CD(0)', CD(0/180*np.pi, ctx))
print('CD(TO)', CD(9.7/180*np.pi, ctx))
print('CD0(0)', CD0(0/180*np.pi, ctx))
print('D_in(0)', 0.5*ctx.rho*w.S*ctx.V**2*CDin(0/180*np.pi, ctx))
print('D_zl(0)', 0.5*ctx.rho*w.S*ctx.V**2*(CDv(0/180*np.pi, ctx) + CDif(0/180*np.pi, ctx) + CDmisc(0/180*np.pi, ctx)))
# clw, clt = CL(0/180*np.pi)
# print('Lw', 0.5*ctx.rho*w.S*ctx.V**2*clw)
# print('Lt', 0.5*ctx.rho*w.S*ctx.V**2*clt)
print('RE', w.Re_mac(ctx))

# working: fx = 4
fx = 3.8
figAR = 1.33144073811

# Airfoil polars
alphad = np.linspace(-10, 18, 100)
alpha = alphad/180*np.pi
Res = np.array([1e6, 1.5e6, 2e6, 2.5e6])
plt.figure(figsize=(fx, fx/figAR))
for Re in Res:
    plt.plot(alphad, w.foil.Cl(alpha, Re), label='Re = %.2E' % Re)
plt.xlabel('Angle of attack [\\SI{}{\\degree}]')
plt.ylabel('Lift coefficient')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('airfoil_cl.pdf')
plt.show()

plt.figure(figsize=(fx, fx/figAR))
for Re in Res:
    plt.plot(w.foil.Cd(alpha, Re), w.foil.Cl(alpha, Re), label='Re = %.2E' % Re)
plt.xlabel('Drag coefficient')
plt.ylabel('Lift coefficient')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('airfoil_polar.pdf')
plt.show()

plt.figure(figsize=(fx, fx/figAR))
for Re in Res:
    plt.plot(alphad, w.foil.Cl(alpha, Re)/w.foil.Cd(alpha, Re), label='Re = %.2E' % Re)
plt.xlabel('Angle of attack [\\SI{}{\\degree}]')
plt.ylabel('Lift-to-drag ratio')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('airfoil_ltd.pdf')
plt.show()

plt.figure(figsize=(fx, fx/figAR))
for Re in Res:
    plt.plot(alphad, w.foil.Cm(alpha, Re), label='Re = %.2E' % Re)
plt.xlabel('Angle of attack [\\SI{}{\\degree}]')
plt.ylabel('Moment coefficient')
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('airfoil_cm.pdf')
plt.show()

# Stall speed
m = np.linspace(580, 700, 100)
Vs = np.zeros_like(m)
CLs_max = np.zeros_like(m)
for i, mm in enumerate(m):
    a, b = CL_stall(mm, ctx)
    CLs_max[i] = a
    Vs[i] = b
plt.figure(figsize=(6.5, 3/0.9*0.94))
plt.plot(m*2.20462262, Vs*1.94384449)
plt.axvline(1322.77357, linestyle='--', c='#8c8b82')
plt.text(1322.77357, 0.5, 'Nominal', color='#8c8b82', ha='right', va='top', rotation=90,
            transform=plt.gca().get_xaxis_transform())
plt.axvline(1543.23584, linestyle='--', c='#8c8b82')
plt.text(1543.23584, 0.5, 'Max.', color='#8c8b82', ha='right', va='top', rotation=90,
            transform=plt.gca().get_xaxis_transform())
plt.xlabel('Takeoff weight [lb]')
plt.ylabel('Stall speed [kt]')
plt.grid()
plt.tight_layout()
plt.savefig('m_vstall.pdf')
plt.show()

plt.plot(m*2.20462262, CLs_max)
plt.show()

# plt.figure(figsize=(7, 3/0.9*0.94))
# alpha = np.linspace(-1, 8, 100)/180*np.pi
# alphad = alpha*180/np.pi
# cdv = CDv(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdv, label='$C_{D,\\mathrm{visc.}}$')
# plt.fill_between(alpha/np.pi*180, cdv, alpha=0.5)
# cdin = cdv + CDin(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdin, label='$+ C_{D,\\mathrm{ind.}}$')
# plt.fill_between(alpha/np.pi*180, cdin, cdv, alpha=0.5)
# cdif = cdin + CDif(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdif, label='$+ C_{D,\\mathrm{int.}}$')
# plt.fill_between(alpha/np.pi*180, cdif, cdin, alpha=0.5)
# cdmisc = cdif + CDmisc(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdmisc, label='$+ C_{D,\\mathrm{misc.}}$')
# plt.fill_between(alpha/np.pi*180, cdmisc, cdif, alpha=0.5)
# plt.axvline(0, linestyle='--', c='#8c8b82')
# plt.text(0, 0.52, 'Cruise', color='#8c8b82', ha='right', va='top', rotation=90,
#             transform=plt.gca().get_xaxis_transform())
# plt.xlabel('Angle of attack [°]')
# plt.ylabel('Drag coefficient $C_D$')
# plt.xticks(np.arange(min(alphad), max(alphad)+1, 1.0))
# plt.yticks(np.arange(0.0, max(cdmisc), 0.005))
# plt.legend(loc='upper left')
# plt.tight_layout()
# plt.savefig('drag.pdf')
# plt.show()

# plt.figure(figsize=(7, 3/0.9*0.94))
# alpha = np.linspace(-1, 8, 100)/180*np.pi
# alphad = alpha*180/np.pi
# cdtot = CD(alpha, ctx)
# cdv = CDv(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdv/cdtot*100, label='$C_{D,\\mathrm{visc.}}$')
# cdin = CDin(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdin/cdtot*100, label='$C_{D,\\mathrm{ind.}}$')
# cdif = CDif(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdif/cdtot*100, label='$C_{D,\\mathrm{int.}}$')
# cdmisc = CDmisc(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdmisc/cdtot*100, label='$C_{D,\\mathrm{misc.}}$')
# plt.axvline(0, linestyle='--', c='#8c8b82')
# plt.text(0, 0.52, 'Cruise', color='#8c8b82', ha='right', va='top', rotation=90,
#             transform=plt.gca().get_xaxis_transform())
# plt.xlabel('Angle of attack [°]')
# plt.ylabel('Percentage of total drag [\\%]')
# plt.xticks(np.arange(min(alphad), max(alphad)+1, 1.0))
# plt.legend(loc='upper left')
# plt.tight_layout()
# plt.savefig('drag_repartition.pdf')
# plt.show()

# plt.figure(figsize=(7, 3/0.9*0.94))
# alpha = np.linspace(-1, 8, 100)/180*np.pi
# alphad = alpha*180/np.pi
# cdin = CDin_test(alpha, ctx)
# plt.plot(alpha/np.pi*180, cdin, label='$C_{D,\\mathrm{ind.}}$')
# plt.axvline(0, linestyle='--', c='#8c8b82')
# plt.text(0, 0.52, 'Cruise', color='#8c8b82', ha='right', va='top', rotation=90,
#             transform=plt.gca().get_xaxis_transform())
# plt.xlabel('Angle of attack [°]')
# plt.ylabel('Percentage of total drag [\\%]')
# plt.xticks(np.arange(min(alphad), max(alphad)+1, 1.0))
# plt.legend(loc='upper left')
# plt.tight_layout()
# plt.savefig('drag_repartition.pdf')
# plt.show()

# Polar
#plt.figure(figsize=(5, 4))
# alpha = np.linspace(-10, 15, 100)/180*np.pi
# alphad = alpha*180/np.pi
# clt = CL_VLM(alpha, ctx)
# cdt = CD(alpha, ctx)
# plt.plot(cdt, clt)
# sampler = len(alphad)//20
# alphads = alphad[::sampler]
# cdts = cdt[::sampler]
# clts = clt[::sampler]
# plt.scatter(cdts, clts, marker='^', color='black')
# for i, val in enumerate(alphads):
#     plt.annotate('\\SI{%.1f}{\\degree}' % (val), (cdts[i], clts[i]), xytext=(-0.6, -0.4), textcoords='offset fontsize', horizontalalignment='right')
# plt.xlabel('Drag coefficient $C_D$')
# plt.ylabel('Lift coefficient $C_L$')
# plt.xlim(xmin=0, xmax=0.05)
# plt.ylim(ymin=-0.6, ymax=1.5)
# plt.grid()
# plt.text(0.002, 0.575 + 0.0275, 'Cruise $C_L$', c='#8c8b82')
# plt.axhline(0.575, linestyle='--', c='#8c8b82')
# plt.savefig('polar.pdf', bbox_inches='tight')
# plt.show()

# LtD
alpha = np.linspace(-4, 10, 100)/180*np.pi
alphad = alpha*180/np.pi
clt = CL_VLM(alpha, ctx)
cdt = CD(alpha, ctx)
plt.plot(alphad, clt/cdt)
sampler = len(alphad)//15
alphads = alphad[::sampler]
cdts = cdt[::sampler]
clts = clt[::sampler]
ltds = clts/cdts
plt.scatter(alphads, ltds, marker='^', color='black')
for i, vclts in enumerate(clts):
    val = np.sqrt(600*9.81 / (0.5 * ctx.rho * 18 * vclts))
    plt.annotate('$%.1f$' % (val * 1.944), (alphads[i], ltds[i]), xytext=(-0.6, -0.4), textcoords='offset fontsize', horizontalalignment='right')
plt.xlabel('Angle of attack [\\SI{}{\degree}]')
plt.ylabel('Lift-to-drag ratio')
plt.axhline(30, linestyle='-', c='#7db928')
plt.axhline(25, linestyle='-', c='#e62d31')
plt.axvline(0, linestyle='--', c='#8c8b82')
plt.text(0, 0.52, 'Cruise', color='#8c8b82', ha='right', va='top', rotation=90,
            transform=plt.gca().get_xaxis_transform())
plt.grid()
plt.savefig('ltd.pdf', bbox_inches='tight')
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