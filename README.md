# Challenges and Opportunities for Single-Sample Network Modeling

This repository contains the code used to run the analyses and create the figures in our paper "Challenges and Opportunities for Single-Sample Network Modeling". A pre-print of the paper is available here: https://www.biorxiv.org/content/10.64898/2026.02.27.708608

The main script is "RunAllAnalyses.m", which calls functions in both the "AnalysisCode" and "NetworkCode" directories. Figures are printed to the "Figures" directory. Intermediate data files are printed to the "Data" directory. All directories are fully populated with the output from running "RunAllAnalyses.m".

The original input data files can be obtained from this zenodo repository: https://doi.org/10.5281/zenodo.18509155.
Prior to executing the "RunAllAnalyses.m" code, these files should be added into the "Data" directory and the "GTExData.zip" file should be un-zipped.

The "NetworkCode" directory contains matlab functions which run the different single-sample network algorithms in various ways.

The "AnalysisCode" directory contains matlab functions which process the data and/or analyze the single-sample networks estimated by different approaches.

The "matlab", "python", and "R" directories contain stand-alone implementations of the four 'linear' methods evaluated in the paper, including LIONESS::PCC, SSN, SWEET, and BONOBO. Insights into scalability and run-times of these implementations are also contained within those directories.
