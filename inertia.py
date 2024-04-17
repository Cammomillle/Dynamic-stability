import numpy as np
from scipy.optimize import fmin
import pos_ellipse as pe
from data import *

rho = 1800

x = np.linspace(0, l_fus, 1000)
dx = np.diff(x)
cx = (x[:-1] + x[1:])/2

def volume():
    a = pe.a_i(cx*1000)/1000
    b = pe.b_i(cx*1000)/1000
    dS = np.pi * a*b
    dV = dS * dx

    return np.sum(dV)

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

def Ac():
    amax = np.max(pe.a) / 1000
    bmax = np.max(pe.b) / 1000
    return np.pi*amax*bmax

def derivative(func):
    """Make a derivative function of a function passed as argument"""
    eps = 0.000001
    return lambda x: (func(x+eps) - func(x-eps))/(2*eps)

def compute_x0():
    # Find x1 (maximum of derivative of body shape; here the lower body line is used)
    dsdx = lambda x: -derivative(pe.ll_i)(x)
    x1 = fmin(dsdx, 3000)[0]/1000

    # Get x0 (Corr. slide 49)
    x0 = l_fus*(0.378 + 0.527*(x1/l_fus))
    
    return x0

print('t_90', thickness(90))
print('amax', amax())
print(volume())
print('Ac', Ac())
print('x0', compute_x0())