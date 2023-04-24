import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

font = {'family': 'times new roman',
        'weight': 'normal',
        'size': 12}

mpl.rc('font', **font)
mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [4, 6]

k = 7.39 * 10**22

plt.figure(1)

v_sw = np.sqrt(1000) # km/s
B_sw = 10**(-8) # teslas
B_s = np.linspace(1, 20, 1000)
a = 5

L_r = k * (B_s)**(2/3) * a * v_sw * B_sw**2

plt.subplot(2,1,1)
plt.plot(B_s, L_r)
plt.xlabel(r"$B_s$ $(G)$")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $B_s$")

v_sw = np.sqrt(1000) # km/s
B_sw = 10**(-8) # teslas
B_s = np.sqrt(20)
a = np.linspace(0.1, 10, 1000)

L_r = k * (B_s)**(2/3) * a * v_sw * B_sw**2

plt.subplot(2,1,2)
plt.plot(a, L_r)
plt.xlabel(r"$a$ $(AU)$")
plt.ylabel(r"$L_r$ (Watts)")
plt.title(r"Scaling of $L_r$ by $a$")

plt.show()
plt.close(1)

