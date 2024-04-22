g = 9.80665
L = 0.0065
T0 = 288.15
M = 0.0289652
R = 8.31446
p0 = 101325
mu0 = 1.789e-5
S = 109.91832

# Calculate pressure with altitude
def pressure(h):
    return p0*(1 - L*h/T0)**(g*M/R/L - 1)

# Calculate the temperature
def temperature(h):
    return T0 - L*h

# Calculate the density 
def density(h):
    return pressure(h)*M/R/T0

def viscosity(h):
    T = temperature(h)
    return mu0*(T/T0)**(3/2) * (T0 + S)/(T + S)