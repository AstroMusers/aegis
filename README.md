# Radio Signatures of Exoplanets

This repository is a comprehensive one containing anything done as 
part of the Radio Signatures of Exoplanets research project run by 
Asaf Kaya and Tansu Daylan. The codes use scientific knowledge acquired
from various research papers and books, most notably the paper 
"Detecting Magnetospheric Radio Emission from Giant Exoplanets" by 
Ashtari et al. 2022.
#
## Python Scripts 
### The Main Prediction Script
The main script that predicts exoplanetary emission quantities is the file
"exoplanet_predictions.py". The scripts implements a Monte-Carlo error propagation
method to be used on exoplanets selected from NASA Exoplanet Archive.

### The Scaling Script
The Python script "scaling.py""
serves the purpose of visualizing the scaling of radio luminosity of
exoplanets with some of their parameters. These are expected not to 
have any direct contribution to any research paper that might be published.

### The Module
The script "radio_module.py" provides a new module with an exoplanet class
and various functions that might be used throughout the project. 

### Expectation Script and the GUI.
The scripts "synthetic_predictions.py" and "GUI_for_predictions.py" deal with
a generated sample of exoplanets. The former is the main script that generates
a random sample of exoplanets by randomly assigning values drawn from
specific probability distributions to various characteristics of exoplanets.
It then plots the expectations of CMI radiation frequency and radio
brightnesses of the drawn sample.
#

## Supplementary Files
The csv files are data files used in the Python scripts. The file
"dataAshtari22.txt" is the result data file from Ashtari 2022, kept here
for comparison. 
#

## Output Files
The txt files "freq.txt", "intens.txt", and "names.txt" are all output tables
that are essentially various representations of the same result obtained in 
the main prediction script.They contain  expected flux density and maximum 
emission frequency data for the exoplanets
chosen from NASA Exoplanet Archive. Every one of these files are sorted 
with respect to different columns, in accordance with the file names.
#

## Results
Resulting plots are gathered in the folder "Results".
# 
## Known Problems
Known problems of the project are summarized in the markdown file 
"known_problems.md".