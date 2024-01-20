## Estimation code for several Commodity pricing models

The models we discuss here possess an affine (or quasi) formula for the future commodity pricing. The estimation method is based on Linear Kalman filtering (and so likelihood maximation). 

The methodology is briefly: we observe future prices (with different maturities) from which we can extrapolate spot prices and other state variables which might help explaining the behavior of the spot prices. Note that in some markets the spot price is not observed. 

The code is developed in Matlab language and consists of the following functions:
a) Starter function: provides randomized starter points which serve as inputs
b) Likelihood-prediction decompotion: function that casts the model in state-space form and then performs filtering.

# References:
1) E. Schwartz (1997) - "The stochastic behavior of commodity prices: Implications for valuation and hedging"
   Contains three models
3) 
