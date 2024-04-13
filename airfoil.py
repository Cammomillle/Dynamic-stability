import numpy as np
from scipy.interpolate import interp1d, RegularGridInterpolator, interp2d
from scipy.optimize import fmin

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