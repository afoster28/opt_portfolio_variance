# opt_portfolio_variance
## Option portfolio variance in CPP

Authors: Adam Foster, Sercan Akbay

The program feeds in a portfolio of options, samples underlying asset prices at maturity, calculates payoffs, portfolio payoffs and summary statistics per portfolio and persists final data and reference data.

Prerequisites: C++ compiler to run `opt_v12.cpp`, input data saved in the same directory (`sample_small.csv`, `sample_medium.csv`, `sample_large.csv`)

The program performs the following processing steps:
* Read csv input containing options
* Parameters given: $S = 5000, μ = 5\\%, σ = 20\\%
* Extract days to maturity of the options
* Calculate
  * t = Days To Maturity/365 = 60/365$
  * $μ_{log} = log(S) + (μ - σ^2/2)t$
  * $σ_{log} = \sqrt{σ^2t}$
  * The variables were converted as above using the Black-Scholes framework, resulting in parameters for the normal distribution of log prices
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
  * $Payout_p = \sum_{q=0}^n Payout_q$
* Calculate mean, variance and standard deviation for every portfolio
* Store summary statistics per portfolio in a results matrix
* Persist the results matrix and reference data (containing prices $S$ and option payouts $Payout_q$)

Program features:
* The workflow is easy to follow with a clear progression from prices $\rightarrow$ option payouts $\rightarrow$ portfolio payouts $\rightarrow$ summary statistics
* Obtaining the underlying price distribution is efficient as it involves sampling from a distribution with known parameters, drawing on Black-Scholes theory and no path dependence for European options
  * This avoids tedious Monte Carlo simulations for price path, reducing the number of parameters and computation time
  * The price distribution is generated once at the start of the process
* Persisting summary statistics `result_matrix.csv` and intermediate price-payout matrix `combined_table.csv` to CSV
* Flexibility: the user can choose a 'complete' (entire input CSV is processed in one go) or 'sequential' run (in stages)
  * The stages in 'Sequential' are portfolios with an incremental number of options, i.e. the first one contains option 1, the second one contains options 1 & 2, etc.
  * Summary statistics are presented on-screen before and after each portfolio
  * Both 'complete' and 'sequential' make use of the same fundamental processing functions
* User-friendliness: in the 'sequential' run the user can freely set a variance cutoff which stops the program if variance exceeds this cutoff, choose whether to proceed to the next portfolio (adding the next option), reset the process
  * The summary statistics and intermediate price-payout matrices reflect these user choices

Potential improvements:
* Model objects could be reallocated to more classes and function scope set more accurately
* Even more efficient approaches could be considered, e.g. sampling based on the outcomes of the normal distribution given its known shape or the covariance approach
* A manual option entry feature for the user
* UI for remote launching of the program
* Integrated and ongoing testing

Testing:
* Active printing during program run: options processed, number of options loaded, portfolio statistics, data saving
* Confirmed shape of the price sample is lognormal
* Spot-checked option payouts as being in the right direction as per option type type and trade directionality
* Reproduced `result_matrix.csv` mean, variance and standard deviation on payout data persisted in `combined_table.csv` under the 'complete' run
* Reconciled `result_matrix.csv`with the expected tables in `output_small.csv`, `output_medium.csv` and `output_large.csv`, showing low percentage differences