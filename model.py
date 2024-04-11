import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator, interp2d
from scipy.optimize import fmin
from data import *

class Airfoil:
    def __init__(self, data):
        interps = {}
    
        angle_min = 1000
        angle_max = -1000
        n = 0

        for k, v in data.items():
            polar = np.genfromtxt(v, delimiter=',')
            polar[:,0] *= np.pi/180
            interps[k] = {
                'cl': interp1d(polar[:,0], polar[:,1], kind='quadratic', bounds_error=False, fill_value='extrapolate'),
                'cd': interp1d(polar[:,0], polar[:,2], kind='quadratic', bounds_error=False, fill_value='extrapolate'),
                'cm': interp1d(polar[:,0], polar[:,4], kind='quadratic', bounds_error=False, fill_value='extrapolate')
            }
            angle_min = min(angle_min, polar[0,0])
            angle_max = max(angle_max, polar[-1,0])
            n = max(n, polar.shape[0])
        
        ns = np.linspace(angle_min, angle_max, n)
        Res = np.array(list(data.keys()))
        Cl_data = np.full((n, len(Res)), np.nan)
        Cd_data = np.full((n, len(Res)), np.nan)
        Cm_data = np.full((n, len(Res)), np.nan)
        for j, Re in enumerate(Res):
            Cl_data[:, j] = interps[Re]['cl'](ns)
            Cd_data[:, j] = interps[Re]['cd'](ns)
            Cm_data[:, j] = interps[Re]['cm'](ns)
        
        Cl = RegularGridInterpolator((ns, Res), Cl_data, bounds_error=False, fill_value=np.nan)
        Cd = RegularGridInterpolator((ns, Res), Cd_data, bounds_error=False, fill_value=np.nan)
        Cm = RegularGridInterpolator((ns, Res), Cm_data, bounds_error=False, fill_value=np.nan)
        self.Cl = lambda a, b: Cl((a, b))
        self.Cd = lambda a, b: Cd((a, b))
        self.Cm = lambda a, b: Cm((a, b))
    
    def alpha0(self, Re):
        return fmin(lambda x: self.Cl(x, Re)**2, 0, disp=False)[0]
    
    def Cl_max(self, Re):
        r = fmin(lambda x: -self.Cl(x, Re), 0, disp=False)[0]
        return self.Cl(r, Re)
    
    def Cl_alpha(self, alpha, Re):
        eps = 1e-6
        a = self.Cl(alpha-eps, Re)
        b = self.Cl(alpha+eps, Re)
        return (b - a)/(2*eps)

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

    def CL(self, ctx):
        return self.mCL(self, ctx)

    def CDv(self, ctx):
        return self.mCDv(self, ctx)
    
    def CDif(self, ctx):
        return self.mCDv(self, ctx)
    
    def CDin(self, ctx):
        return self.mCDin(self, ctx)

    def CD(self, ctx):
        return self.mCDv(self, ctx) + self.mCDif(self, ctx) + self.mCDin(self, ctx)

def WingCL_alpha_DATCOM(w, ctx):
    M = ctx.V / 330
    beta = np.sqrt(1 - M**2)
    k = w.foil.Cl_alpha(ctx.alpha_e, w.Re_mac(ctx)) / (2*np.pi)
    return (2*np.pi*w.AR)/(2 + np.sqrt((w.AR*beta/k)**2 * (1 + np.tan(w.sweep(0.5))**2 / beta**2) + 4))

def WingCL_DATCOM(w, ctx):
    CL_alpha = w.CL_alpha(ctx)
    alpha0 = w.foil.alpha0(w.Re_mac(ctx))

    return CL_alpha * (ctx.alpha - alpha0)

def WingCDif_Zero(w, ctx):
    return 0

class WingCDin_Oswald:
    def __init__(self, e):
        self.e = e

    def __call__(self, w, ctx):
        CL = w.CL(ctx)
        return CL**2 / (np.pi * self.e * w.AR)

class WingCDv_Strip:
    def __init__(self, n, half=True):
        self.n = n
        self.half = half

    def __call__(self, w, ctx):
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
        Cds = w.foil.Cd(ctx.alpha, Re)*dS

        # Calculate wingtip drag
        S_wl = (w.c_tip + w.wl_tip) * w.wl_h / 2
        Cd_wl = 2 * w.foil.Cd(ctx.beta, Re[-1]) * S_wl / w.S

        return np.sum(Cds)/(w.S/d) + Cd_wl

class Fuselage:
    def CDin(self, alpha):
        
    



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


class Ctx:
    pass

ctx = Ctx()
ctx.alpha = 0/180 * np.pi
ctx.alpha_e = 0
ctx.beta = 0
ctx.V = 32
ctx.rho = 1.22
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
h.mCL_alpha = WingCL_alpha_DATCOM
h.mCL = WingCDif_Zero
h.mCDv = WingCDv_Strip(50, False)
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
w.mCL_alpha = WingCL_alpha_DATCOM
w.mCL = WingCL_DATCOM
w.mCDv = WingCDv_Strip(50)
w.mCDif = WingCDif_Zero
w.mCDin = WingCDin_Oswald(0.9335)

print(w.CD(ctx) + v.CD(ctx)*v.S/w.S + h.CD(ctx)*h.S/w.S)
def polar():
    alpha = np.linspace(-8, 8, 30)/180 * np.pi
    CL = np.zeros_like(alpha)
    CD = np.zeros_like(alpha)
    for i, a in enumerate(alpha):
        ctx.alpha = a
        CL[i] = w.CL(ctx)
        CD[i] = w.CD(ctx)

    plt.plot(CD, CL)
    plt.show()