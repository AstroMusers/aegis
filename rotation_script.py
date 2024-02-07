import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

names, types, per, flag = np.genfromtxt("rotation.csv", delimiter=",", dtype=str, skip_header=1, unpack=True)

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
    if flag[i] == "1":
        p = time_string_to_days(per[i])
    else:
        p = float(per[i])
    periods.append(p)

# print(periods)

log_periods = np.log10(periods)

# print(log_periods)

kde = gaussian_kde(log_periods)

num_samples = 500
samples = kde.resample(size=num_samples)

true_samples = 10 ** samples

true_samples = true_samples.T
samples = samples.T

plt.figure(1)
plt.hist(log_periods, bins=29, color="teal", edgecolor="black")
plt.xlabel('$\log_{10}{\mathrm{(Rotational\,\, Period / days)}}$')
plt.ylabel("Occurence")
plt.title("Histogram of Solar System Bodies' Rotation Periods")
# Create histogram
plt.figure(2)
plt.hist(samples, bins=30, color='skyblue', edgecolor='black')
plt.xlabel('$\log_{10}{\mathrm{(Simulated\,\, Rotational\,\, Period / days)}}$')
plt.ylabel('Density (arbitrary units)')
plt.title('Histogram of Random Samples')
#
# # Show plot
plt.show()



