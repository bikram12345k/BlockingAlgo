# BlockingAlgo
Implementation of an approximation algorithm for blocking of an experimental design

This repository includes the code to implement the algorithms for blocking. 
The general algorithm is implemented in 'bottlenecknbpv2.R'

Specific applications of the algorithm and its modification can be found in the
folders.  The folders are created based on the presentation of the results in
the paper and its supplement.  

Section 3 compares the algorithm to 'nbpMatching' 

Section 4 compares 8 algorithms including that of Moore (2012), Higgins et al. (2016) 
and Karmakar (2018).

Section 5 shows the calculation of the balance values in the NHANES data.

Supplement folder contains 
(1) Code to replicate the results showing the efficincy of the proposed algo over 
standard algo in treatment effect estimation.
(2) More simulation results complementing Section 4 results (k=3 and k=8).
(3) Analysis of the the NHANES data.

# 9/10/2020
Additional codes are included that adds outcome model to the simulation cases
and compares different methods for treatment effect estimation.

# 5/13/2021
The organization of the results have been changed. The earlier code has not been removed.
The code to reproduce the results after this reorganization is given in the folder
"REORGANIZATION". The file names correspond to the related table.
