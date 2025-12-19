# Challenges and Opportunities for Single-Sample Network Modeling

This repository contains the code used to run the analyses and create the figures in our paper "Challenges and Opportunities for Single-Sample Network Modeling".

The main script is "RunAllAnalyses.m", which calls functions in both the "AnalysisCode" and "NetworkCode" directories. Figures are printed to the "Figures" directory and intermediate files are printed to the "Data" directory.

The "NetworkCode" directory contains matlab functions which run the different single-sample network algorithms in various ways.

The "AnalysisCode" directory contains matlab functions which process the data and/or analyze the single-sample networks estimated by different approaches.

The "python" directory contains two scripts to complement the "GenerateSSCorr.m" function in "NetworkCode". The python script "GenerateSSCorr.py" simultaneously runs four different single-sample network approaches based on Pearson correlation (LIONESS::PCC, SSN, SWEET, BONOBO) on the input data. The python script "run_SSN_method.py" allows the user to select one of these single-sample network approaches to run on the input data.
