# Equity Securities
This should update now.

```@docs
VLQuantitativeFinancePackage.EulerMaruyamaMethod
VLQuantitativeFinancePackage.MyOrnsteinUhlenbeckModel
VLQuantitativeFinancePackage.MyHestonModel
```

## State Space Models of Return
Update this section with some groooovy text.

```@docs
VLQuantitativeFinancePackage.MySisoLegSHippoModel
VLQuantitativeFinancePackage.build(model::Type{MySisoLegSHippoModel}, data::NamedTuple)
VLQuantitativeFinancePackage.solve(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64})
VLQuantitativeFinancePackage.prediction(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64,1};
    S::Int64 = 10, B::Float64 = 40.0, α::Float64 = 0.25, β::Float64 = 0.10)
```