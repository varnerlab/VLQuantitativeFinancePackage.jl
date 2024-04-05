using VLQuantitativeFinancePackage

# parameters -
r̄ = 0.049; # risk free rate
Δt = (1/365); # step-size, assuming 365-trading days per year
T = 62*Δt; # duration of the contract (units: years)
IV = 0.4005; # implied volatility - we are estimating this
Sₒ = 79.50; # share price at contract purchase (units: USD/share)
K = 75.0; # strike price for the MU contract (units: USD/share)
h = 248; # number of levels in the binomial tree
𝒫 = 2.93; # contract premium

# build contract
test_american_put_contract_model = build(MyAmericanPutContractModel, (
        K = K, sense = 1, copy = 1, DTE = T, IV = 0.1, premium = 𝒫));

# compute the implied volatility -
test_IV = estimate_implied_volatility(test_american_put_contract_model, Sₒ = Sₒ, h = h, r̄ = r̄)