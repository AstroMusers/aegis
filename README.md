# Aegis

This is a comprehensive repository containing everything done as 
part of the research project run by Asaf Kaya and Dr. Tansu Daylan
on the Cyclotron Maser Instability (CMI) driven radio emissions
from exoplanets. The codes use scientific knowledge acquired
from various research papers and books.

The most essential prior file that needs to exist in the same directory with the
python scripts is a .csv file containing data from NASA Exoplanet Archive with
the required parameters. The name of this file is "NASA[sample].csv", where
[sample] is usually the date the data was downloaded.

## Python Scripts 

### The Modules
The script "radio_module.py" provides a new module with an exoplanet class
and various functions that are used throughout the project. The script
"rotation_script.py" provides an additional module to randomly sample spin angular
momenta of exoplanets from the solar system distribution. Both of these modules
are essential to be able to run the main prediction script.

### The Wind Script
"wind.py" imports data from a .csv file obtained from the NASA Exoplanet
Archive with the required parameters, whose name is specified in the script.
It calculates the stellar wind temperatures and speed at the orbit
of the exoplanets in the sample assuming isothermal Parker wind (Parker 1958).
The results are written out to be used in the main prediction script into file
"wind_info[sample].date"

### The Main Prediction Script
The main script that predicts exoplanetary emission characteristics is the file
"exoplanet_predictions.py". The script implements a Monte-Carlo error propagation
method to be used on exoplanets selected from NASA Exoplanet Archive, in order to 
quantify uncertainty in the results.
The script requires the existence of the exoplanet archive table in .csv format.


### Extra Visual Routines
The programs "parker_spiral.py", "sketch.py" and "visibility.py" are used to create
three of the figures
in the paper: the perpendicular component of the IMF in the orbit of an exoplanet
(currently tau Boo b), the schematic drawing of a magnetized exoplanet, and the visibility
figure containing all-sky maximum elevation and time-spent-above-20-degrees maps
for the considered telescopes, respectively. The first one requires the existence
of wind data from an exoplanet, currently tau Boo b, to be present in the directory
in .npz format. This file is currently "taub_wind.npz". The second is an independent script,
while the latter requires the existence of the file "obs_table.csv" in the same 
working directory. This file is created within "extract.py", which extracts the
necessary information for the opportune targets determined from
"exoplanet_predictions.py"

### Obsolete Scripts
#### The Scaling Script
The Python script "scaling.py""
serves the purpose of visualizing the scaling of radio luminosity of
exoplanets with some of their parameters. These are expected not to 
have any direct contribution to any research paper that might be published.
#### Expectation Script and the GUI
The scripts "synthetic_predictions.py" and "GUI_for_predictions.py" deal with
a generated sample of exoplanets. The former is the main script that generates
a random sample of exoplanets by randomly assigning values drawn from
specific probability distributions to various characteristics of exoplanets.
It then plots the expectations of CMI radiation frequency and radio
brightnesses of the drawn sample.

## Output Files
The resulting radio flux densities and maximum emission frequencies
from "exoplanet_predictions.py" are written
out in two manners. First, the potentially visible candidates determined from
the magnetic and kinetic RBLs, and the candidates determined from an integration
of both RBLs are separated, sorted by their names, frequencies, and flux densities
separately into nine .txt files in "old_result_tables". 
Secondly, only the "both" RBL methods results for all exoplanets are 
written out to csv files with the uncertainties included in "Output Tables".
In this case, the results for all exoplanets can be found in "all.csv", while
there also exist four different files that are subsets of this file divided from
the expected frequencies.

## Known Problems
Known problems of the project are summarized in the markdown file 
"known_problems.md".