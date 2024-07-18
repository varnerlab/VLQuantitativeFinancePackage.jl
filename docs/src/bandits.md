# Bandit Problems
Fill me in

## Episilon Greedy Sampling
```@docs
VLQuantitativeFinancePackage.MyEpsilonSamplingBanditModel
VLQuantitativeFinancePackage.build(model::Type{MyEpsilonSamplingBanditModel}, data::NamedTuple)
```

## Ticker picker problem example
Fill me in

```@docs
VLQuantitativeFinancePackage.MyTickerPickerWorldModel
VLQuantitativeFinancePackage.build(model::Type{MyTickerPickerWorldModel}, data::NamedTuple)
VLQuantitativeFinancePackage.sample(model::MyEpsilonSamplingBanditModel, world::AbstractWorldModel)
VLQuantitativeFinancePackage.preference
```

Risk aware bandit problem example
```@docs
VLQuantitativeFinancePackage.MyTickerPickerRiskAwareWorldModel
VLQuantitativeFinancePackage.build(model::Type{MyTickerPickerRiskAwareWorldModel}, data::NamedTuple)
```