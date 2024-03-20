# Calcium-Imaging
Calcium imaging files from ImageJ to MatLab to further analysis


A simple MatLab pipeline for analysing a time series obtained from ImageJ

Created by Valtteri Vuolanto and Fikret Emre Kapucu, NeuroGroup, Faculty of Medicine and Health Technology, Tampere University, Tampere, Finland
and implemented first in https://www.researchsquare.com/article/rs-3399985/v1.
Cite the article and github URL https://github.com/VVuolanto/Calcium-Imaging for further usage.

Input is a folder containing .csv files obtained from ImageJ image processing software, 
with use of Time Series Analyzer v3.0 -plugin.
The data is formated as frames in rows and each column being a ROI, i.e. a detected neuron from the imaged region.

Each neuron is processed separately, smoothing data to dispose of small background variations. The data is normalized in regards to baseline fluoresence with dF/F0.
The data is also corrected by subtracting a linear approximation; this may not be desireable for all datasets but those with significant photo-bleaching.
Resulting signal is analyzed for peaks, and multiple parameters (derived) from the peaks are stored in a .xlsx file in desired location.
