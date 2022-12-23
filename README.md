This repo contains an example R code used in the Monte Carlo experiments of the paper "Robust function-on-function interaction regression"
# Authors
Ufuk Beyaztas, Han Lin Shang and Abhijit Mandal
# Procedure
The main file is run. The rest files include the auxilary functions.
# Packages
library('MASS')
library('fda')
library('fda.usc')
library('rrcov')
library('ROCR')
library('data.table')
library('regRSM')
library('FRegSigCom')
library('mvtnorm')
# Note
Please note that the package FRegSigCom is not on CRAN anymore. Install the FRegSigCom package using the provided .tar.gz file to perform LQ method
