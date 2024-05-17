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
We model [United States Treasury securities](https://www.treasurydirect.gov), e.g., [Treasury bills](https://www.treasurydirect.gov/marketable-securities/treasury-bills/), [Treasury notes](https://www.treasurydirect.gov/marketable-securities/treasury-notes/), and [Treasury bonds](https://www.treasurydirect.gov/marketable-securities/treasury-bonds/) using the `MyUSTreasuryZeroCouponBondModel` and `MyUSTreasuryCouponSecurityModel` types, which are subtypes of the `AbstractTreasuryDebtSecurity` abstract type.  

```@docs
VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel
VLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel
```

We construct `MyUSTreasuryZeroCouponBondModel` and `MyUSTreasuryCouponSecurityModel` objects by specifying the maturity date, the face value of the bond the effective annual interest rate and other data that are specific to the bill, note or bond in a `build` function where the data is 
specified in a [NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) format. 

For example, let's compute the price and cashflow for a `T = 20-yr` bond, with a coupon rate of `c = 1.750%`, a yield (discount rate) `rate = 1.850%`, two coupon payments per year, i.e., $\lambda = 2$ and a face (par) value of $V_{P}$ = `100 USD`

```julia
test_bond = build(MyUSTreasuryCouponSecurityModel, (
    T = 20.0, rate = 0.01850, coupon = 0.01750, λ = 2, par = 100.0
)) |> discount_model;
```


```@docs
VLQuantitativeFinancePackage.price
VLQuantitativeFinancePackage.strip
```
