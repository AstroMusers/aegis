import tkinter as tk
from tkinter import ttk
import subprocess


# Function to save the most recent parameters to a configuration file
def save_recent_parameters():
    recent_parameters = [
        entry1.get(), entry2.get(), entry3.get(), entry4.get(), entry5.get(),
        entry6.get(), entry7.get(), entry8.get(), entry9.get(), entry10.get(),
        entry11.get()
    ]

    with open("recent_parameters.txt", "w") as file:
        file.write("\n".join(recent_parameters))


# Function to load the most recent parameters from the configuration file
def load_recent_parameters():
    try:
        with open("recent_parameters.txt", "r") as file:
            recent_parameters = file.read().splitlines()

        # Populate the entry fields with the most recent parameters
        entry1.delete(0, tk.END)
        entry1.insert(tk.END, recent_parameters[0])
        entry2.delete(0, tk.END)
        entry2.insert(tk.END, recent_parameters[1])
        entry3.delete(0, tk.END)
        entry3.insert(tk.END, recent_parameters[2])
        entry4.delete(0, tk.END)
        entry4.insert(tk.END, recent_parameters[3])
        entry5.delete(0, tk.END)
        entry5.insert(tk.END, recent_parameters[4])
        entry6.delete(0, tk.END)
        entry6.insert(tk.END, recent_parameters[5])
        entry7.delete(0, tk.END)
        entry7.insert(tk.END, recent_parameters[6])
        entry8.delete(0, tk.END)
        entry8.insert(tk.END, recent_parameters[7])
        entry9.delete(0, tk.END)
        entry9.insert(tk.END, recent_parameters[8])
        entry10.delete(0, tk.END)
        entry10.insert(tk.END, recent_parameters[9])
        entry11.delete(0, tk.END)
        entry11.insert(tk.END, recent_parameters[10])
    except FileNotFoundError:
        # Handle the case when the configuration file does not exist
        pass


def run_predictions():
    # Get the input values
    save_recent_parameters()
    a_min = entry1.get()
    a_max = entry2.get()
    B_mean = entry3.get()
    B_sd = entry4.get()
    M_mean = entry5.get()
    M_sd = entry6.get()
    Mdot_mean = entry7.get()
    Mdot_sd = entry8.get()
    D_mean = entry9.get()
    D_sd = entry10.get()
    sample_size = entry11.get()
    emission_type = burst.get()

    # Validate and convert the input values to float
    try:
        a_min = float(a_min)
        a_max = float(a_max)
        B_mean = float(B_mean)
        B_sd = float(B_sd)
        M_mean = float(M_mean)
        M_sd = float(M_sd)
        Mdot_mean = float(Mdot_mean)
        Mdot_sd = float(Mdot_sd)
        D_mean = float(D_mean)
        D_sd = float(D_sd)
        sample_size = int(sample_size)
        emission_type = int(emission_type)
    except ValueError:
        # Handle any invalid input
        print("Invalid input. Please enter numeric values.")
        return

    # Execute the synthetic_predictions.py script with the input values as arguments
    subprocess.call(['python', 'synthetic_predictions.py', str(a_min), str(a_max), str(B_mean), str(B_sd),
                     str(M_mean), str(M_sd), str(Mdot_mean), str(Mdot_sd), str(D_mean), str(D_sd), str(sample_size), str(emission_type)])


# Create the main window
window = tk.Tk()
window.title("Exoplanet Predictions")
window.geometry("500x400")

burst = tk.IntVar()

# Default values for the variables
default_a_min = -2
default_a_max = 2
default_B_mean = 10
default_B_sd = 3
default_M_mean = 3
default_M_sd = 0.8
default_Mdot_mean = 100
default_Mdot_sd = 30
default_D_mean = 1000
default_D_sd = 300
default_size = 150

# Create entry fields for input variables with default values
style = ttk.Style()
style.theme_use('classic')  # Choose theme
style.configure("TEntry", padding=5, font=("Arial", 12))

label1 = ttk.Label(window, text="Minimum of semi-major axis (logAU):")
label1.grid(row=0, column=0, sticky="w")
entry1 = ttk.Entry(window, style="Custom.TEntry")
entry1.insert(tk.END, default_a_min)
entry1.grid(row=0, column=1)

label2 = ttk.Label(window, text="Maximum of semi-major axis (logAU)::")
label2.grid(row=1, column=0, sticky="w")
entry2 = ttk.Entry(window, style="Custom.TEntry")
entry2.insert(tk.END, default_a_max)
entry2.grid(row=1, column=1)

label3 = ttk.Label(window, text="Mean of Magnetic Field Strength (Gauss):")
label3.grid(row=2, column=0, sticky="w")
entry3 = ttk.Entry(window, style="Custom.TEntry")
entry3.insert(tk.END, default_B_mean)
entry3.grid(row=2, column=1)

label4 = ttk.Label(window, text="SD of Magnetic Field Strength (Gauss):")
label4.grid(row=3, column=0, sticky="w")
entry4 = ttk.Entry(window, style="Custom.TEntry")
entry4.insert(tk.END, default_B_sd)
entry4.grid(row=3, column=1)

label5 = ttk.Label(window, text="Mean of the host star's mass (Solar Masses):")
label5.grid(row=4, column=0, sticky="w")
entry5 = ttk.Entry(window, style="Custom.TEntry")
entry5.insert(tk.END, default_M_mean)
entry5.grid(row=4, column=1)

label6 = ttk.Label(window, text="SD of the host star's mass (Solar Masses):")
label6.grid(row=5, column=0, sticky="w")
entry6 = ttk.Entry(window, style="Custom.TEntry")
entry6.insert(tk.END, default_M_sd)
entry6.grid(row=5, column=1)

label7 = ttk.Label(window, text="Mean of the host star's mass loss rate (10^-15 Solar Masses per year):")
label7.grid(row=6, column=0, sticky="w")
entry7 = ttk.Entry(window, style="Custom.TEntry")
entry7.insert(tk.END, default_Mdot_mean)
entry7.grid(row=6, column=1)

label8 = ttk.Label(window, text="SD of the host star's mass loss rate (10^-15 Solar Masses per year):")
label8.grid(row=7, column=0, sticky="w")
entry8 = ttk.Entry(window, style="Custom.TEntry")
entry8.insert(tk.END, default_Mdot_sd)
entry8.grid(row=7, column=1)

label9 = ttk.Label(window, text="Mean of the host star distance (light years):")
label9.grid(row=8, column=0, sticky="w")
entry9 = ttk.Entry(window, style="Custom.TEntry")
entry9.insert(tk.END, default_D_mean)
entry9.grid(row=8, column=1)

label10 = ttk.Label(window, text="SD of the host star distance (light years):")
label10.grid(row=9, column=0, sticky="w")
entry10 = ttk.Entry(window, style="Custom.TEntry")
entry10.insert(tk.END, default_D_sd)
entry10.grid(row=9, column=1)

# Create an empty label as a separator
empty_label = ttk.Label(window, text="")
empty_label.grid(row=10, column=0, columnspan=2)

label11 = ttk.Label(window, text="Sample size:")
label11.grid(row=11, column=0, sticky="e")
entry11 = ttk.Entry(window, style="Custom.TEntry")
entry11.insert(tk.END, default_size)
entry11.grid(row=11, column=1)

checkbox = ttk.Checkbutton(window,
                           text="Burst?",
                           variable=burst,
                           onvalue=1,
                           offvalue=0)
checkbox.grid(row=12, column=0)

load_recent_parameters()

# Define the custom entry style with a grayish background color and a black frame
style.configure("Custom.TEntry", background="#f0f0f0", bordercolor="black", relief="solid", borderwidth=1)

# Create a button to run the predictions
button = ttk.Button(window, text="Run Predictions", command=run_predictions)
button.grid(row=13, column=0, columnspan=3)

# Start the GUI event loop
window.mainloop()
