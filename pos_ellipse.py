import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Données fournies
pos = np.array([0, 1, 10, 49.7895, 100, 200, 500, 800, 1200, 1500, 2000, 2500, 3002.5, 3501.5, 4008.3, 4523.5, 5000, 5500.9, 6000.7, 6500.7, 7002.3, 7501.8, 8001.8, 8291.3, 8348.8, 8398.8195822750000000, 8456.7, 8483.5, 8493.3, 8496.6, 8500])

init_y = -13.23019572106
mid_line_y = - init_y + np.array([-13.23019572106, -13.21975540270, -13.21975540270, -13.21975540270, -17.20931390022, -19.49432025010, 3.02332490648, 35.63977669865, 82.77825439860, 119.34105574323, 179.66236119902, 237.13061229928, 293.78568043059, 358.93525040696, 421.16602551569, 450.37260033156, 472.93008397190, 492.85579427838, 508.36555842451, 519.08643982637, 530.39311213552, 544.20315203958, 558.49332050205, 566.76999067656, 573.10530063734, 588.21316145546, 629.91433799045, 656.68758648633, 666.51462115367, 669.84328774590, 672.56544475503])

a = np.array([0, 41.2293178110000000, 101.7208773390000000, 202.0197467410000000, 262.5784903300000000, 357.1550745440000000, 559.5276297800000000, 709.2603349740000000, 853.8425682250000000, 925.6424095370000000, 984.1090572820000000, 972.7901193840000000, 889.4498456540000000, 721.2439047670000000, 545.5412154440000000, 439.7599213890000000, 361.3731824340000000, 297.9164244360000000, 254.5001170470000000, 231.9328491770000000, 219.5167001120000000, 212.3522730070000000, 206.0134717040000000, 202.3393156740000000, 192.2269678250000000, 164.2344037860000000, 83.4082459830000000, 31.0506434330000000, 11.8329539750000000, 5.3234337550000000, 0])/2
b = np.array([0, 15.7797338350000000, 48.9897948530000000, 99.8417782760000000, 117.5776996860000000, 150.6586809480000000, 233.8712523590000000, 294.0952844290000000, 340.2693935840000000, 350.0151878220000000, 342.6148862910000000, 320.4911108820000000, 283.3643634100000000, 232.3063443860000000, 156.6160699430000000, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100])

ml_i = interp1d(pos, mid_line_y, kind='quadratic')
a_i = interp1d(pos, a, kind='quadratic')
b_i = interp1d(pos, b, kind='quadratic')
ll_i = lambda x: ml_i(x) - a_i(x)
hl_i = lambda x: ml_i(x) + a_i(x)
bh_i = lambda x: 2*a_i(x)
bw_i = lambda x: 2*b_i(x)

# Calcul des positions verticales des points supérieurs et inférieurs
y_top = mid_line_y + a
y_low = mid_line_y - a

# Ajustement des valeurs à x=0
y_top -= y_top[0]
y_low -= y_low[0]

# Interpolation des points supérieurs
f_top = interp1d(pos, y_top, kind='linear')

# Interpolation des points inférieurs
f_low = interp1d(pos, y_low, kind='linear')

# Interpolation des points supérieurs pour la vue du dessus
f_b_top = interp1d(pos, -b, kind='linear')

# Interpolation des points inférieurs pour la vue du dessus
f_b_low = interp1d(pos, b, kind='linear')

# Création des nouvelles positions pour un tracé plus précis
new_pos = np.linspace(pos[0], pos[-1], 1000)

# Calcul des valeurs interpolées
new_y_top = f_top(new_pos)
new_y_low = f_low(new_pos)
new_b_top = f_b_top(new_pos)
new_b_low = f_b_low(new_pos)

def xflr_text(n_frame, n_lines):
    text = """ 
# This file defines a body geometry
# The frames are defined from nose to tail
# The numer of sidelines is defined by the number of points of the first frame
# Each of the next frames should have the same number of points as the first
# For each frame, the points are defined for the right half of the body, 
# in the clockwise direction aft looking forward

Body Name

BODYTYPE
 1        # Flat Panels (1) or NURBS (2)

OFFSET
0.0     0.0     0.0     #Total body offset (Y-coord is ignored)

"""
    frames = np.linspace(pos[0], pos[-1], n_frame)
    lines = np.linspace(np.pi/2, -np.pi/2, n_lines)

    for frame in frames:
        x = frame
        a = a_i(x)
        b = b_i(x)
        ml = ml_i(x)

        r = a*b/np.sqrt((b*np.sin(lines))**2 + (a*np.cos(lines))**2)
        r[np.isnan(r)] = 0.0
        y = r*np.cos(lines)
        z = ml + r*np.sin(lines)
        plt.plot(y, z, label='%0.2f' % x)

        text += 'FRAME\n'
        for y, z in zip(y, z):
            x_sp = ' '*(4 if x < 0 else 5)
            y_sp = ' '*(9 if y < 0 else 10)
            z_sp = ' '*(7 if z < 0 else 8)

            text += '%s%.7f%s%.7f%s%.7f\n' % (x_sp, x/1000, y_sp, y/1000, z_sp, z/1000)
        text += '\n'

    return text

def side_area():
    A = 0
    for i in range(len(new_y_top)):
        A += (new_y_top[i]-new_y_low[i])
    return A

def max_height_width_ratio():
    max_height = np.max(new_y_top-new_y_low)
    max_width = np.max(new_b_low-new_b_top)
    return max_height/max_width

#m_lin = W_fus/(2*np.pi*g)   # assumption of a symmetry of revolution for the mass computation


#print("Side area: ", side_area()/10**6, "m^2")
#print("h/w: ", max_height_width_ratio())


