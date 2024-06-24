# Derivative Securities
Derivatives are financial instruments based on the value of other assets like commodities, stocks, or market indexes. Options are a type of derivative that uses stock as its underlying asset. They are contractual agreements giving the buyer the right, but not the obligation, to execute a transaction at a later date. 

## Options contracts
Options contracts are financial instruments that give the holder the right to buy or sell an asset at a predetermined price on or before a specific date. The two main types of options are call options and put options. 

* Call options give the holder the right to buy an asset at a predetermined price on or before the expiration date.
* Put options give the holder the right to sell an asset at a predetermined price on or before the expiration date.

### European option contracts
European options are contracts that give the holder the right to buy or sell an asset at a predetermined price on the expiration date. 
```@docs
VLQuantitativeFinancePackage.MyEuropeanCallContractModel
VLQuantitativeFinancePackage.MyEuropeanPutContractModel
VLQuantitativeFinancePackage.build(model::Type{MyEuropeanCallContractModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyEuropeanPutContractModel}, data::NamedTuple)
```

### American option contracts
American options are contracts that give the holder the right to buy or sell an asset at a predetermined price at any time before the expiration date.
```@docs
VLQuantitativeFinancePackage.MyAmericanCallContractModel
VLQuantitativeFinancePackage.MyAmericanPutContractModel
VLQuantitativeFinancePackage.build(model::Type{MyAmericanCallContractModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyAmericanPutContractModel}, data::NamedTuple)
```

## Contracts at expiration
```@docs
VLQuantitativeFinancePackage.payoff
VLQuantitativeFinancePackage.profit
```

## European contract premiums
```@docs
VLQuantitativeFinancePackage.MyBlackScholesContractPricingModel
VLQuantitativeFinancePackage.premium(contract::MyEuropeanCallContractModel, 
    model::MyBlackScholesContractPricingModel; sigdigits::Int64 = 4)
VLQuantitativeFinancePackage.premium(contract::MyEuropeanPutContractModel, 
    model::MyBlackScholesContractPricingModel; sigdigits::Int64 = 4)
```

## American contract premiums
```@docs
VLQuantitativeFinancePackage.MyAdjacencyBasedCRREquityPriceTree
VLQuantitativeFinancePackage.MyCRRLatticeNodeModel
VLQuantitativeFinancePackage.build(model::Type{MyAdjacencyBasedCRREquityPriceTree}, data::NamedTuple)
VLQuantitativeFinancePackage.populate(model::MyAdjacencyBasedCRREquityPriceTree; 
    Sₒ::Float64 = 100.0, h::Int64 = 1)
VLQuantitativeFinancePackage.premium
```

## Implied volatility
```@docs
VLQuantitativeFinancePackage.estimate_implied_volatility
```

## The Greeks
[The Greeks](https://en.wikipedia.org/wiki/en:Greeks_(finance)?variant=zh-tw) quantify the sensitivity of an option's premium to various factors. `Delta`, `theta`, `vega`, `rho` and `gamma` are the most widely used Greeks, measuring an option's sensitivity to changes in the underlying asset's share price, the time decay, implied volatility, the risk-free rate and the rate of change in delta, respectively. 

```@docs
VLQuantitativeFinancePackage.delta
VLQuantitativeFinancePackage.theta
VLQuantitativeFinancePackage.vega
VLQuantitativeFinancePackage.rho
VLQuantitativeFinancePackage.gamma
VLQuantitativeFinancePackage.vwap
```