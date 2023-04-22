import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

v_sw = np.linspace(100,1000,1000) # km/s
R_mp = 10000 #km
B_sw = 6 * 10**(-9)
mu_0 = 4 * np.pi * 10**(-7)

L_r = 3 * 10**(-3) * np.pi * R_mp**2 * v_sw * B_sw**2 / (2*mu_0)

plt.figure(1)
plt.plot(v_sw, L_r)
plt.show()