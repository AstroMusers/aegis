import matplotlib.pyplot as plt
import numpy as np

# Generate sample data
x = np.linspace(0, 10, 100)
y = np.sin(x)

# Create a figure and axis
fig, ax = plt.subplots()

# Plot your data
im = ax.imshow([y], aspect='auto', cmap='viridis', extent=[0, 10, 0, 1])

# Add colorbar
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('Colorbar Label')

# Adjust layout to remove whitespace
plt.tight_layout()

# Show the plot
plt.show()
