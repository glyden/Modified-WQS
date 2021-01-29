# Penalized weighted sum regression for chemical mixtures
This repository contains the code to fit a modified form of Weighted Quantile Sum regression that does not require data to be split into training and test sets for weight estimation and inference. This method was introduced in "A permutation-based approach to inference for weighted sum regression with correlated chemical mixtures" by Lyden, et al.

The contents of this repository are as follows:

functions.R includes functions that are helpful in generating correlated data.

method-functions.R includes functions to fit modified WQS, i.e., estimate mixture component weights using penalized constrained regressions and run a permutation test for the overall mixture effect.

example.R applies method-functions.R to estimate weights and overall mixture effects in a simulated dataset.

Simulations/ contains code to run all simulation studies described in Lyden, et al.
