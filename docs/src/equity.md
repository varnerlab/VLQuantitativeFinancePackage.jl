# Equity Securities
This should update now.

## Computing returns
Fill me in.

```@docs
VLQuantitativeFinancePackage.log_growth_matrix
```

## Lattice Models
Fill me in

## Geometric Brownian Motion
Fill me in

## Advanced Stochastic Differential Equation Models
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