This directory contains a script that to complements the "GenerateSSCorr.m" function in "main/NetworkCode".

The R script "GenerateSSCorr.R" contains five functions:

GenerateSSCorr() simultaneously runs four different single-sample network approaches based on Pearson correlation (LIONESS::PCC, SSN, SWEET, BONOBO) on the input data.

The functions GenerateLIONESS(), GenerateSSN(), GenerateSWEET(), and GenerateBONOBO() allows the user to run only one of these single-sample network approaches on the input data.

The "scalability_tests" directory contains an analysis showing how these methods scale based on the dimensions of the input data. All tests were run using R version 4.5.2 on Windows 11 Pro, with an Intel(R) Core(TM) Ultra 7 165U processor and 32 GB of installed RAM.
