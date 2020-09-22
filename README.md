# Penalized weighted sum regression for chemical mixtures
This repository contains the code to fit a modified form of Weighted Quantile Sum regression that does not require an analyst to split data into training and test sets for weight estimation and inference, as introduced in "A permutation-based approach to inference for weighted sum regression with correlated chemical mixtures" by Lyden, et al.

The contents are as follows:

functions.R includes functions that are helpful in generating correlated data.

method-functions.R has all the functions to fit our method, i.e., run the penalized constrained regressions that estimate the weights and execute a permutation test for the overall mixture effect.

example.R shows how to call the above scripts to estimate weights and overall mixture effects in a simulated dataset.
