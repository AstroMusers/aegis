import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["figure.autolayout"] = True
plt.rcParams['figure.figsize'] = [5, 3]

rc = {"font.family": "times new roman",
      "font.size": 11,
      "mathtext.fontset": "stix"}
plt.rcParams.update(rc)

mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["navy", "darkred", "b"])

k = 7.39 * 10 ** 22

v_sw = 400  # km/s
B_sw = 10 ** (-8)  # teslas
B_s = np.linspace(1, 50, 1000)
a = np.linspace(0.1, 10, 1000)

B_s, a = np.meshgrid(B_s, a)

L_r = k * B_s ** (2 / 3) * a ** (2 / 3) * v_sw * B_sw ** 2

fig, ax = plt.subplots(1, 1)
# cp = ax.contourf(X, Y, L_r)
cp = plt.imshow(L_r, cmap=plt.cm.viridis, origin="lower",
                extent=[B_s.min(), B_s.max(), a.min(), a.max()], aspect="auto")
fig.colorbar(cp, label=r"$L_r$ (Watts)")  # Add a colorbar to a plot
ax.set_title('Scaling of $L_r$ by $B_s$ and $a$')
ax.set_xlabel(r'$B_s$ (G)')
ax.set_ylabel(r'$a$ (AU)')
plt.show()
