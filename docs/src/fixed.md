# Fixed Income Securities
Fixed income securities are financial instruments that pay a fixed amount of interest over a specified period of time. The most common fixed income securities are bonds, which are issued by governments, municipalities, and corporations. The fixed income market is one of the largest financial markets in the world, and it plays a critical role in the global economy.

## Discounting models
In the `VLQuantitativeFinancePackage` we allow computing the present value of a cash flow stream using a discounting model. 
The present value of a cash flow stream can be computed using a `ContinuousDiscountingModel` or a `DiscreteDiscountingModel`.

```@docs
VLQuantitativeFinancePackage.DiscreteCompoundingModel
VLQuantitativeFinancePackage.ContinuousCompoundingModel
```

## Pricing models
We model [United States Treasury securities](), e.g., [Treasury bills](https://www.treasurydirect.gov/marketable-securities/treasury-bills/), [Treasury notes](https://www.treasurydirect.gov/marketable-securities/treasury-notes/), and [Treasury bonds](https://www.treasurydirect.gov/marketable-securities/treasury-bonds/) using the `MyUSTreasuryZeroCouponBondModel` and `MyUSTreasuryCouponSecurityModel` types, which are subtypes of the `AbstractTreasuryDebtSecurity` abstract type.  

```@docs
VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel
VLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel
```




```@docs
VLQuantitativeFinancePackage.price
VLQuantitativeFinancePackage.strip
```
