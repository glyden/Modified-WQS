# Penalized weighted sum regression for chemical mixtures
This repository contains the code to fit a modified form of Weighted Quantile Sum regression that does not require data to be split into training and test sets for weight estimation and inference. This method was introduced in "A permutation-based approach to inference for weighted sum regression with correlated chemical mixtures" by Lyden, et al.

The contents are as follows:

functions.R includes functions that are helpful in generating correlated data.

method-functions.R has all the functions to fit our method, i.e., fit the penalized constrained regressions to estimate the weights and run a permutation test for the overall mixture effect.

example.R shows how to use the functions to estimate weights and overall mixture effects in a simulated dataset.
