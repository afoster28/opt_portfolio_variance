# opt_portfolio_variance
##Option portfolio variance in CPP

Authors: Adam Foster, Sercan Akbay

The program feeds in a portfolio of options, samples underlying asset prices at maturity, calculates payoffs, portfolio payoffs and summary statistics per portfolio and persists final data and reference data.

Prerequisites: C++ compiler, input data in the same directory as `opt_v12.cpp` (e.g. sample_small.csv)

The program performs the following processing steps:
* Read csv input containing options
* Calculate $μ_{log} = log(S) + (μ - σ^2/2)t$