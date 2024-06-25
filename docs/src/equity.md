# Equity Securities
This should update now.

## Computing returns
Fill me in.

```@docs
VLQuantitativeFinancePackage.log_growth_matrix
VLQuantitativeFinancePackage.vwap
```

## Lattice Models
```@docs
VLQuantitativeFinancePackage.sample(model::MyBinomialEquityPriceTree, L::Int64; 
    number_of_paths::Int64 = 100)
```

## Geometric Brownian Motion
```@docs
VLQuantitativeFinancePackage.MyGeometricBrownianMotionEquityModel
VLQuantitativeFinancePackage.MyMultipleAssetGeometricBrownianMotionEquityModel
VLQuantitativeFinancePackage.build(model::Type{MyGeometricBrownianMotionEquityModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyMultipleAssetGeometricBrownianMotionEquityModel}, data::NamedTuple)
VLQuantitativeFinancePackage.sample(model::MyGeometricBrownianMotionEquityModel, data::NamedTuple; 
    number_of_paths::Int64 = 100)
VLQuantitativeFinancePackage.sample(model::MyMultipleAssetGeometricBrownianMotionEquityModel, data::NamedTuple; 
    number_of_paths::Int64 = 100)
VLQuantitativeFinancePackage.sample_endpoint(model::MyGeometricBrownianMotionEquityModel, data::NamedTuple; 
    number_of_paths::Int64 = 100)
```

## Advanced Stochastic Pricing and Return Models
```@docs
VLQuantitativeFinancePackage.EulerMaruyamaMethod
VLQuantitativeFinancePackage.MyOrnsteinUhlenbeckModel
VLQuantitativeFinancePackage.MyHestonModel
VLQuantitativeFinancePackage.build(model::Type{MyOrnsteinUhlenbeckModel}, data::NamedTuple)
VLQuantitativeFinancePackage.build(model::Type{MyHestonModel}, data::NamedTuple)
VLQuantitativeFinancePackage.solve(model::AbstractAssetModel, tspan::NamedTuple,
    initialconditions::AbstractArray; method::AbstractStochasticSolverModel = EulerMaruyamaMethod(), 
    N::Int64 = 100)
```

## Structured State-Space Models
Update this section with some groooovy text.

```@docs
VLQuantitativeFinancePackage.MySisoLegSHippoModel
VLQuantitativeFinancePackage.build(model::Type{MySisoLegSHippoModel}, data::NamedTuple)
VLQuantitativeFinancePackage.solve(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64})
VLQuantitativeFinancePackage.prediction(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64,1};
    S::Int64 = 10, B::Float64 = 40.0, α::Float64 = 0.25, β::Float64 = 0.10)
VLQuantitativeFinancePackage.estimate_hippo_parameters(model::MySisoLegSHippoModel, tspan::NamedTuple, 
    signal::Array{Float64}; method = LBFGS())
```