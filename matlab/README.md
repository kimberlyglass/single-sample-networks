This directory contains scripts that complement the "GenerateSSCorr.m" function in "main/NetworkCode".

Within "method_code" directory there are five functions written in MATLAB:

GenerateSSCorr() simultaneously runs four different single-sample network approaches based on Pearson correlation (LIONESS::PCC, SSN, SWEET, BONOBO) on the input data.

The functions GenerateLIONESS(), GenerateSSN(), GenerateSWEET(), and GenerateBONOBO() allows the user to run only one of these single-sample network approaches on the input data.

The "scalability_tests" directory contains an analysis showing how these methods scale based on the dimensions of the input data. All tests were run using MATLAB R2026a on Windows 11 Pro, with an Intel(R) Core(TM) Ultra 7 165U processor and 32 GB of installed RAM.
