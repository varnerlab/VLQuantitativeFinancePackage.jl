# Fixed Income Securities
Fixed income securities are financial instruments that pay a fixed amount of interest over a specified period of time. The most common fixed income securities are bonds, which are issued by governments, municipalities, and corporations. The fixed income market is one of the largest financial markets in the world, and it plays a critical role in the global economy.


## Discounting models
In the `VLQuantitativeFinancePackage` we allow computing the present value of a cash flow stream using a discounting model. The present value of a cash flow stream can be computed using a [`ContinuousDiscountingModel`](@ref) or a [`DiscreteDiscountingModel`](@ref).

```@docs
VLQuantitativeFinancePackage.DiscreteCompoundingModel
VLQuantitativeFinancePackage.ContinuousCompoundingModel
```

## Treasury models
We model [United States Treasury securities](https://www.treasurydirect.gov), e.g., [Treasury bills](https://www.treasurydirect.gov/marketable-securities/treasury-bills/), [Treasury notes](https://www.treasurydirect.gov/marketable-securities/treasury-notes/), and [Treasury bonds](https://www.treasurydirect.gov/marketable-securities/treasury-bonds/) using 
the [`MyUSTreasuryZeroCouponBondModel`](@ref) and [`MyUSTreasuryCouponSecurityModel`](@ref) types, which are subtypes of the `AbstractTreasuryDebtSecurity` abstract type.  

```@docs
VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel
VLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel
```

## Computing the price of a fixed income Treasury security
We construct [`MyUSTreasuryZeroCouponBondModel`](@ref) and [`MyUSTreasuryCouponSecurityModel`](@ref) objects by specifying the maturity date, the face value of the bond the effective annual interest rate and other data that are specific to the bill, note or bond in a `build` function where the data is specified in a [NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) format. 


```@docs
VLQuantitativeFinancePackage.price
```

## Separating the principal and interest payments
[Registered Interest and Principal of Securities (STRIPS) bonds](https://en.wikipedia.org/wiki/United_States_Treasury_security#STRIPS) are a unique type of fixed-income investment instrument that provides investors with an alternative way to access the coupon payments of Treasury securities. STRIPS bonds are created by separating a Treasury securities coupon and principal components and trading them as individual  zero-coupon securities. This process allows investors to purchase and trade the coupon or principal components separately, providing greater flexibility in managing their investment portfolios.


```@docs
VLQuantitativeFinancePackage.strip
```