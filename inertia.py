import numpy as np
from scipy.optimize import fmin
import pos_ellipse as pe
from data import *

rho = 1800

x = np.linspace(0, l_fus, 1000)
dx = np.diff(x)
cx = (x[:-1] + x[1:])/2

def amax():
    a = pe.a_i(x*1000)/1000
    b = pe.b_i(x*1000)/1000
    S = np.pi*a*b
    return np.max(S)

def mass(t):
    a = pe.a_i(cx*1000)/1000
    b = pe.b_i(cx*1000)/1000
    dS = np.pi * (a*b - (a-t)*(b-t))
    dV = dS * dx

    V = np.sum(dV)
    M = V*rho
    print('center of gravity', np.sum(cx * dV / V))
    return M

def thickness(m):
    return fmin(lambda t: (mass(t) - m)**2, 0.005, disp=False)[0]

print('t_90', thickness(90))
print('amax', amax())
