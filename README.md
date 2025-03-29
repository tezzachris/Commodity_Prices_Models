## Commodity Stochastic Factor Models 

The models we consider here all possess a quasi-affine formula for the pricing of commodity futures. 

Which factors? Typically in the literature: log-spot commodity price, convenience yield, interest rate, volatility of log-price, long-term mean of log-price
Which commodities? We considered mainly crude-oil, copper, natural gas

# The Models
- One-two-three factor in [1]
- Three factor in [3]
- Four factor in [2]
- Four factor in [4]

## The Estimation Method

In commodity markets future prices at several maturities are quoted. From this prices we can extrapolate the factors/latent variables I described above. The estimation method is usually based on the Kalman filter, and so on likelihood maximation. 
Remark: in some markets the spot commodity price is not quoted. 

In general terms, the code is developed in Matlab language and consists of the following functions:
(a) Starter: provides randomized starter points for the coefficients 
(b) Likelihood prediction error decomposition: casts the model into its state-space form, performs Kalman filtering and computes the log-likelihood.
(c) Obtain the MLEs by repeating (a) and (b) to avoid local minima

Note: We rely on Matlab built-in optimizers (see "fminsearch" and "fmincon")

# References:
[1] E. Schwartz (1997) - "The stochastic behavior of commodity prices: Implications for valuation and hedging" \
[2] X.S. Yan (2002) - "Valuation of commodity derivatives in a new multi-factor model"\
[3] W.K. Hughen (2010) - "A maximal affine stochastic volatility model of oil prices"\
[4] S. Spinler & Schone (2017) - "A four-factor stochastic volatility model of commodity prices"\
[5] Ballestra, Tezza (2025) - "A multi-factor model for improved commodity pricing: Calibration and an application to the oil market"
