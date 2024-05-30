# Derivative Securities
Fill me in

## European and American Options Contracts
```@docs
VLQuantitativeFinancePackage.MyEuropeanCallContractModel
VLQuantitativeFinancePackage.MyEuropeanPutContractModel
VLQuantitativeFinancePackage.MyAmericanCallContractModel
VLQuantitativeFinancePackage.MyAmericanPutContractModel
VLQuantitativeFinancePackage.build(model::Type{MyEuropeanCallContractModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyEuropeanPutContractModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyAmericanCallContractModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyAmericanPutContractModel}, data::NamedTuple)
```

## Payoff and profit of contracts at expiration
```@docs
VLQuantitativeFinancePackage.payoff
VLQuantitativeFinancePackage.profit
```

## Computing European contract premiums
```@docs
VLQuantitativeFinancePackage.MyBlackScholesContractPricingModel
VLQuantitativeFinancePackage.premium(contract::MyEuropeanCallContractModel, 
    model::MyBlackScholesContractPricingModel; sigdigits::Int64 = 4)
VLQuantitativeFinancePackage.premium(contract::MyEuropeanPutContractModel, 
    model::MyBlackScholesContractPricingModel; sigdigits::Int64 = 4)
```

## Computing American contract premiums
```@docs
VLQuantitativeFinancePackage.MyAdjacencyBasedCRREquityPriceTree
VLQuantitativeFinancePackage.MyCRRLatticeNodeModel
VLQuantitativeFinancePackage.build(model::Type{MyAdjacencyBasedCRREquityPriceTree}, data::NamedTuple)
VLQuantitativeFinancePackage.populate(model::MyAdjacencyBasedCRREquityPriceTree; 
    Sₒ::Float64 = 100.0, h::Int64 = 1)
VLQuantitativeFinancePackage.premium
```

## Estimating the Implied Volatility of American contracts
```@docs
VLQuantitativeFinancePackage.estimate_implied_volatility
```