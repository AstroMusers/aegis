import matplotlib.pyplot as plt
import numpy as np

# # Generate sample data
# x = np.linspace(0, 10, 100)
# y = np.sin(x)
#
# # Create a figure and axis
# fig, ax = plt.subplots()
#
# # Plot your data
# im = ax.imshow([y], aspect='auto', cmap='viridis', extent=[0, 10, 0, 1])
#
# # Add colorbar
# cbar = plt.colorbar(im, ax=ax)
# cbar.set_label('Colorbar Label')
#
# # Adjust layout to remove whitespace
# plt.tight_layout()
#
# # Show the plot
# plt.show()


import matplotlib.pyplot as plt

# Sample plot
x = range(10)
y = [i**2 for i in x]

plt.plot(x, y)

# Annotate with text and arrows
plt.annotate(
    'Observable from ground',
    xy=(5, 50),  # The position where you want the text
    xytext=(5, 60),  # The position of the text
    ha='center',  # Horizontal alignment of the text
    # arrowprops=dict(arrowstyle='<->', lw=1.5),  # Double-headed arrow
    fontsize=12,  # Size of the text
)

# Draw the horizontal line (arrows)
plt.annotate(
    '',
    xy=(1, 50),  # Left end of the arrow
    xytext=(9, 50),  # Right end of the arrow
    arrowprops=dict(arrowstyle='<->', lw=1.5)  # Double-headed arrow
)

plt.show()