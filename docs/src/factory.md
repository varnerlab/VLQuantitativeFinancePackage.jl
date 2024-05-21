# Factory
We provide a factory function `build` to construct most of the types used in this package.

```@docs
VLQuantitativeFinancePackage.build(model::Type{MyUSTreasuryZeroCouponBondModel}, data::NamedTuple)
```

```@docs
VLQuantitativeFinancePackage.build(model::Type{MyUSTreasuryCouponSecurityModel}, data::NamedTuple)
```

```@docs
VLQuantitativeFinancePackage.build(model::Type{MyOrnsteinUhlenbeckModel}, data::NamedTuple)
```

```@docs
VLQuantitativeFinancePackage.build(model::Type{MySisoLegSHippoModel}, data::NamedTuple)
```

```@docs
VLQuantitativeFinancePackage.build(model::Type{MyHestonModel}, data::NamedTuple)
```