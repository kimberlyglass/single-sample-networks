This directory contains a script that to complements the "GenerateSSCorr.m" function in "main/NetworkCode".

The R script "GenerateSSCorr.R" contains five functions:

GenerateSSCorr() simultaneously runs four different single-sample network approaches based on Pearson correlation (LIONESS::PCC, SSN, SWEET, BONOBO) on the input data.

The functions GenerateLIONESS(), GenerateSSN(), GenerateSWEET(), and GenerateBONOBO() allows the user to run only one of these single-sample network approaches on the input data.
