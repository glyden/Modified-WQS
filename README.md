# Penalized weighted sum regression for chemical mixtures
This repository contains the code to fit a modified form of Weighted Quantile Sum regression that allows for the entire dataset to be used in both weight estimation and inference, as introduced in "A permutation-based approach to inference for weighted sum regression with correlated chemical mixtures" by Lyden, et al.

functions.R includes functions that are helpful in generating correlated data.

method-functions.R has all the functions to fit our method, including those based on glmnet to run the penalized constrained regressions that estimate the weights, as well as a function to perform a permutation test for inference on the overall mixture effect.

example.R is a brief example of how to call the above scripts and use them to estimate weights and overall mixture effects in a simulated dataset.
