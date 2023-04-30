import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

rc = {"font.family" : "times new roman",
      "font.size" : 11,
      "mathtext.fontset" : "stix"}
plt.rcParams.update(rc)

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [4, 9]

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["navy", "darkred", "b"])

R_j = 69911 # Jupiter radius in kilometers
mu_0 = 4 * np.pi * 10**(-7) # Vacuum Permeability in SI units


plt.figure(1)

v_sw = np.linspace(1,1000,1000) # km/s
R_mp = 10*R_j
B_sw = 10**(-8)

L_r = 3 * 10**(-3) * np.pi * (R_mp * 10**3)**2 * (v_sw * 10**3) * B_sw**2 / (2*mu_0)

plt.subplot(3,1,1)
plt.plot(v_sw, L_r)
plt.xlabel(r"$v_{sw}$ (km/s)")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $v_{sw}$")

v_sw = 400 # km/s
R_mp = np.linspace(0.1, 100, 1000 ) #Jupiter radius
B_sw = 10**(-8)

L_r = 3 * 10**(-3) * np.pi * (R_mp*R_j* 10**3)**2 * (v_sw * 10**3) * B_sw**2 / (2*mu_0)

plt.subplot(3,1,2)
plt.plot(R_mp, L_r)
plt.xlabel(r"$R_{mp}$ $(R_J)$")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $R_{mp}$")

v_sw = 400 # km/s
R_mp = 10*R_j #km
B_sw = np.linspace(1, 100, 1000)

L_r = 3 * 10**(-3) * np.pi * (R_mp * 10**3)**2 * (v_sw * 10**3) * (B_sw * 10**(-9))**2 / (2*mu_0)

plt.subplot(3,1,3)
plt.plot(B_sw, L_r)
plt.xlabel(r"$B_{sw}$ (nT)")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $B_{sw}$")

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.show()

plt.close(1)
