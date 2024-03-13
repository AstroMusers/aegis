import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from numba import njit
from numba import jit


names, types, per = np.genfromtxt("rotation.csv", delimiter=",", dtype=str, skip_header=1, skip_footer=2, unpack=True, usecols=(0,1,2))
flag, mass, radius = np.genfromtxt("rotation.csv", delimiter=",", dtype=float, skip_header=1, skip_footer=2, unpack=True, usecols=(3,4,5))

def time_string_to_days(time_str):
    # Split the time string into components
    components = time_str.split()

    # Extract hours, minutes, seconds, and milliseconds
    hours = int(components[0][:-1])
    minutes = int(components[1][:-1])
    seconds = int(components[2][:-1])

    # Calculate the total time in seconds
    total_seconds = (hours * 3600) + (minutes * 60) + seconds

    # Convert seconds to days
    days = total_seconds / (24 * 3600)

    # Return days with three decimal places
    return round(days, 3)

# Example usage
time_str = "9h 55m 29s 710ms"
days = time_string_to_days(time_str)
# print(days)

periods = []
for i in range(len(names)):
    if flag[i] == 1:
        p = time_string_to_days(per[i])
    else:
        p = float(per[i])
    periods.append(p)

# print(periods)
periods = np.array(periods)

MJ = mass[np.where(names == "Jupiter")][0]
RJ = radius[np.where(names == "Jupiter")][0]
print(MJ, RJ)

mass /= MJ
radius /= RJ
omega = 2*np.pi/periods

wJ = omega[np.where(names == "Jupiter")][0]
omega /= wJ

momenta = 2 / 5 * mass * radius**2 * omega
momenta /= momenta[np.where(names == "Jupiter")][0]
log_momenta = np.log10(momenta)

# print(log_periods)

kde = gaussian_kde(log_momenta)

num_samples = 500
samples = kde.resample(size=num_samples)

true_samples = 10 ** samples

true_samples = true_samples.T
samples = samples.T

plt.figure(1)
plt.hist(log_momenta, bins=29, color="teal", edgecolor="black")
plt.xlabel('$\log_{10}{\mathrm{(Spin \, Angular \, Momentum \, [L_J])}}$')
plt.ylabel("Occurence")
plt.title("Histogram of Solar System Bodies' Spin Angular Momenta")
# Create histogram
plt.figure(2)
plt.hist(samples, bins=30, color='skyblue', edgecolor='black')
plt.xlabel('$\log_{10}{\mathrm{(Simulated \, Spin \, Angular \, Momentum \, [L_J])}}$')
plt.ylabel('Density (arbitrary units)')
plt.title('Histogram of Random Samples')
#
# # Show plot
plt.show()


def moment_sampler():
    return 10 ** kde.resample(1)[0][0]
