---
title: "readme"
author: "Vishwesha Guttal"
date: "11/22/2016"
output: html_document
---

##Codes and Data for "Alternative stable states and spatial indicators of critical slowing down along a spatial gradient in a savanna ecosystem", 2017, GEB (DOI: http://dx.doi.org/10.1111/geb.12570)

###Codes

1) A C++ code (spacetime.cpp) for the spatially-explicit model to show validity of space-for-time substitution. Written by Sabiha Majumder. 

To run this code, use any latest linux OS. We have run this on Ubuntu 12.04. It is expected to work on other OS as well, if you have the compilers, but we haven’t tested. To run this code, type the following command on the command prompt. 

> g++ -o sptime spacetime.cpp 

This will compile the C++ code and provide an executable named sptime. Now execute this on the command prompt:

> ./sptime &

This will run the code in the backgound. In the default settings, it will output simulation data matrix corresponding to Fig 2 of the main text in a single file. This data is basically a huge matrix of size N x NL (see Box 1 and Appendix S1 for more details). 

2) Spatial Indicator codes to analyse simulation data (by Sabiha Majumder)

a) A C++ code (cg-indicators.cpp) to calculate three spatial indicators (Spatial variance, spatial skewness and spatial correlation). Put this code in the same folder that contains the simulation output of code #1 above. Follow the same procedure as in #1) above to create the executable, but with different names for C++ code and the executable.

b) R code (low-frequency-spectra.R) to calculate spatial low frequency spectrum. This code calls the function powerspectrum from the other file names as powerspectrum.R.  Put these codes in the same folder that contains the simulation output of code #1 above.

3) Spatial indicator (real_sp_indicators.m) codes for transect data (by Amit Agrawal)

This codes takes data from a transect. User specifies certain parameters like dimension of matrix and coarse-graining length. The code outputs the spatial mean, spatial variance, spatial skewness, spatial correlation and spatial DFT at low frequency.  

4) Spatial indicator codes (null_sp_indicators.m) for null data (by Amit Agrawal)

This code first generates null model data and then analyses the spatial indicators as it does for real data in #3. 
5) Spatial Indicators along rainfall gradient (reap_sp_indicators_rainfall.m) code for transect and null model data (by Amit Agrawal).

NOTE: For the analysis corresponding to Fig 3, we have used build-in functions of R/matlab and have explained the procedure in the Methods section. 

###Data

1) Vegetation data is available on the website of Prof. Denne N Reed (http://www.dennereed.org/), and is based on the paper Reed, D. N., et al. "The spatial distribution of vegetation types in the Serengeti ecosystem: the influence of rainfall and topographic relief on vegetation patch characteristics." Journal of Biogeography 36.4 (2009): 770-782.

2) Rainfall data is available from Dr Michael Coughenour (Michael.Coughenour@colostate.edu).

###Contact us
For any queries, contact guttal@ces.iisc.ernet.in

###Reference
Stephanie Eby*, Amit Agrawal*, Sabiha Majumder, Andew Dobson and Vishwesha Guttal, 2017, Alternative stable states and spatial indicators of critical slowing down along a spatial gradient in a savanna ecosystem,  Global Ecology and Biogeography, DOI: http://dx.doi.org/10.1111/geb.12570 (*These authors contributed equally to this work).

###Blogpost
Here is a blogpost on our lab website on this article: https://goo.gl/EyHoSK
