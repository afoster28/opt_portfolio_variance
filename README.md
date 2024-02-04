# opt_portfolio_variance
##Option portfolio variance in CPP

Authors: Adam Foster, Sercan Akbay

The program feeds in a portfolio of options, samples underlying asset prices at maturity, calculates payoffs, portfolio payoffs and summary statistics per portfolio and persists final data and reference data.

Prerequisites: C++ compiler, input data in the same directory as `opt_v12.cpp` (e.g. sample_small.csv)

The program performs the following processing steps:
* Read csv input containing options
* Extract days to maturity of the options
* Calculate
  * $S = 5000, μ = 5\\%, σ = 20\\%, t = Days To Maturity/365 = 60/365$
  * $μ_{log} = log(S) + (μ - σ^2/2)t$
  * $σ_{log} = \sqrt{σ^2t}$
  * The variables were converted as above using the Black-Scholes formula, resulting in parameters for normal distribution of log prices
* Set up random number generator for sampling
* Sample log price 10,000 times
* Exponentiate log prices to get prices
  * The resulting distribution of prices follows a lognormal distibution
* Create a matrix of option payouts: calculate payout $Payout_q$ for every option for every price
  * $Payout C_i = max(0, S - K)$
  * $Payout P_i = max(0, K - S)$
  * $Short \rightarrow Payout_d = - Payout_i$
  * $Payout_q = Payout_d \cdot Quantity$
* Construct a matrix of portfolio payouts: sum option payouts in portfolio $Payout_p$ for every price where portfolios consist of options up to and including the latest option in the iteration
  * $Payout_p = \left( \sum_{q=0}^n Payout_q \right)$
* Calculate mean, variance and standard deviation for every portfolio
* Store summary statistics per portfolio in a results matrix
* Persist the results matrix and reference data (containing prices $S$ and option payouts $Payout_q$)

