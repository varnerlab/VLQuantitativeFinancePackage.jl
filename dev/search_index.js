var documenterSearchIndex = {"docs":
[{"location":"derivatives/#Derivative-Securities","page":"Derivatives","title":"Derivative Securities","text":"","category":"section"},{"location":"derivatives/","page":"Derivatives","title":"Derivatives","text":"Fill me in","category":"page"},{"location":"equity/#Equity-Securities","page":"Equity","title":"Equity Securities","text":"","category":"section"},{"location":"equity/","page":"Equity","title":"Equity","text":"Fill me in","category":"page"},{"location":"fixed/#Fixed-Income-Securities","page":"Fixed Income","title":"Fixed Income Securities","text":"","category":"section"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"Fixed income securities are financial instruments that pay a fixed amount of interest over a specified period of time. The most common fixed income securities are bonds, which are issued by governments, municipalities, and corporations. The fixed income market is one of the largest financial markets in the world, and it plays a critical role in the global economy.","category":"page"},{"location":"fixed/#Discounting-models","page":"Fixed Income","title":"Discounting models","text":"","category":"section"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"In the VLQuantitativeFinancePackage we allow computing the present value of a cash flow stream using a discounting model.  The present value of a cash flow stream can be computed using a ContinuousDiscountingModel or a DiscreteDiscountingModel.","category":"page"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"VLQuantitativeFinancePackage.DiscreteCompoundingModel\nVLQuantitativeFinancePackage.ContinuousCompoundingModel","category":"page"},{"location":"fixed/#VLQuantitativeFinancePackage.DiscreteCompoundingModel","page":"Fixed Income","title":"VLQuantitativeFinancePackage.DiscreteCompoundingModel","text":"struct DiscreteCompoundingModel <: AbstractCompoundingModel\n\nImmutable type that represents discrete compounding. This type has no fields and is passed as an argument to various functions to indicate that discrete compounding should be used in calculations.\n\n\n\n\n\n","category":"type"},{"location":"fixed/#VLQuantitativeFinancePackage.ContinuousCompoundingModel","page":"Fixed Income","title":"VLQuantitativeFinancePackage.ContinuousCompoundingModel","text":"struct ContinuousCompoundingModel <: AbstractCompoundingModel\n\nImmutable type that represents continuous compounding.  This type has no fields and is passed as an argument to various functions to indicate that continuous compounding should be used in calculations.\n\n\n\n\n\n","category":"type"},{"location":"fixed/#Pricing-models","page":"Fixed Income","title":"Pricing models","text":"","category":"section"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"We model United States Treasury securities, e.g., Treasury bills, Treasury notes, and Treasury bonds using the MyUSTreasuryZeroCouponBondModel and MyUSTreasuryCouponSecurityModel types, which are subtypes of the AbstractTreasuryDebtSecurity abstract type.  ","category":"page"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel\nVLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel","category":"page"},{"location":"fixed/#VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel","page":"Fixed Income","title":"VLQuantitativeFinancePackage.MyUSTreasuryZeroCouponBondModel","text":"mutable struct MyUSTreasuryZeroCouponBondModel <: AbstractTreasuryDebtSecurity\n\nA mutable struct that represents a U.S. Treasury zero coupon bond. \n\nFields\n\npar::Float64: Par value of the bond\nrate::Union{Nothing, Float64}: Annual interest rate\nT::Union{Nothing,Float64}: Duration in years, measured as a 365 day or a 52 week year\nprice::Union{Nothing, Float64}: Price of the bond or note\nn::Int: Number of compounding periods per year (typically 2)\ncashflow::Union{Nothing, Dict{Int,Float64}}: Cashflow dictionary where the key is the period and the value is the discounted cashflow in a period\ndiscount::Union{Nothing, Dict{Int,Float64}}: Discount factor dictionary where the key is the period and the value is the discount factor in that period\n\n\n\n\n\n","category":"type"},{"location":"fixed/#VLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel","page":"Fixed Income","title":"VLQuantitativeFinancePackage.MyUSTreasuryCouponSecurityModel","text":"mutable struct MyUSTreasuryCouponSecurityModel <: AbstractTreasuryDebtSecurity\n\nA mutable struct that represents a U.S. Treasury coupon bond.  This type of security (note or bond) pays the holder of the note (or bond) a fixed interest rate at regular intervals over the life of the instrument.\n\nFields\n\npar::Float64: Par value of the bond\nrate::Union{Nothing, Float64}: Annualized effective discount rate\ncoupon::Union{Nothing, Float64}: Coupon interest rate\nT::Union{Nothing,Float64}: Duration in years of the note or bond, measured as a 365 day or a 52 week year\nλ::Int: Number of coupon payments per year (typically 2)\nprice::Union{Nothing, Float64}: Price of the bond or note\ncashflow::Union{Nothing, Dict{Int,Float64}}: Cashflow dictionary where the key is the period and the value is the discounted cashflow in a period\ndiscount::Union{Nothing, Dict{Int,Float64}}: Discount factor dictionary where the key is the period and the value is the discount factor in that period\n\n\n\n\n\n","category":"type"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"We construct MyUSTreasuryZeroCouponBondModel and MyUSTreasuryCouponSecurityModel objects by specifying the maturity date, the face value of the bond the effective annual interest rate and other data that are specific to the bill, note or bond in a build function where the data is  specified in a NamedTuple format. ","category":"page"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"For example, let's compute the price and cashflow for a T = 20-yr bond, with a coupon rate of c = 1.750%, a yield (discount rate) rate = 1.850%, two coupon payments per year, i.e., lambda = 2 and a face (par) value of V_P = 100 USD","category":"page"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"test_bond = build(MyUSTreasuryCouponSecurityModel, (\n    T = 20.0, rate = 0.01850, coupon = 0.01750, λ = 2, par = 100.0\n)) |> discount_model;","category":"page"},{"location":"fixed/","page":"Fixed Income","title":"Fixed Income","text":"VLQuantitativeFinancePackage.price\nVLQuantitativeFinancePackage.strip","category":"page"},{"location":"fixed/#VLQuantitativeFinancePackage.price","page":"Fixed Income","title":"VLQuantitativeFinancePackage.price","text":"price(model::MyUSTreasuryCouponSecurityModel, compounding::T) -> MyUSTreasuryCouponSecurityModel where T <: AbstractCompoundingModel\n\nThe price(...) function computes the price of a MyUSTreasuryCouponSecurityModel instance using a discrete or continuous compounding model.\n\nArguments\n\nmodel::MyUSTreasuryCouponSecurityModel: an instance of the MyUSTreasuryCouponSecurityModel type.\ncompounding::T: an instance of a type that is a subtype of AbstractCompoundingModel, i.e., a discrete or continuous compounding model.\n\nReturns\n\nMyUSTreasuryCouponSecurityModel: an updated instance of the MyUSTreasuryCouponSecurityModel type with the price computed using the compounding model.\n\n\n\n\n\nprice(model::MyUSTreasuryZeroCouponBondModel, compounding::T) -> MyUSTreasuryZeroCouponBondModel where T <: AbstractCompoundingModel\n\nThe price(...) function computes the price of a MyUSTreasuryZeroCouponBondModel instance using a discrete or continuous compounding model.\n\nArguments\n\nmodel::MyUSTreasuryZeroCouponBondModel: an instance of the MyUSTreasuryZeroCouponBondModel type.\ncompounding::T: an instance of a type that is a subtype of AbstractCompoundingModel, i.e., a discrete or continuous compounding model.\n\nReturns\n\nMyUSTreasuryZeroCouponBondModel: an updated instance of the MyUSTreasuryZeroCouponBondModel type with the price computed using the compounding model.\n\n\n\n\n\n","category":"function"},{"location":"fixed/#VLQuantitativeFinancePackage.strip","page":"Fixed Income","title":"VLQuantitativeFinancePackage.strip","text":"strip(model::MyUSTreasuryCouponSecurityModel) -> Dict{Int, MyUSTreasuryZeroCouponBondModel}\n\nStrips the coupon and par value payments from a parent coupon security. \n\nThe strip(...) function takes a model::MyUSTreasuryCouponSecurityModel of the security we wish to strip and returns a Dictionary  holding MyUSTreasuryZeroCouponBondModel instances created from the parent security.  The keys of the dictionary correspond to the temporal index of the created security.\n\n\n\n\n\n","category":"function"},{"location":"#VLQuantitativeFinancePackage.jl","page":"Home","title":"VLQuantitativeFinancePackage.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The VLQuantitativeFinancePackage.jl package is a Julia package that provides a collection of functions and types useful for quantitative finance. The package is designed to be simple and easy to use, and it is suitable for students, researchers, and practitioners in the area of quantitative finance.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed by running the following command in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://github.com/varnerlab/VLQuantitativeFinancePackage.jl.git\")","category":"page"},{"location":"#Disclaimer-and-Risks","page":"Home","title":"Disclaimer and Risks","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This content is offered solely for training and informational purposes. No offer or solicitation to buy or sell securities or derivative products or any investment or trading advice or strategy is made, given, or endorsed by the teaching team. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Trading involves risk. Carefully review your financial situation before investing in securities, futures contracts, options, or commodity interests. Past performance, whether actual or indicated by historical tests of strategies, is no guarantee of future performance or success. Trading is generally inappropriate for someone with limited resources, investment or trading experience, or a low-risk tolerance.  Only risk capital that is not required for living expenses.","category":"page"},{"location":"","page":"Home","title":"Home","text":"You are fully responsible for any investment or trading decisions you make. Such decisions should be based solely on evaluating your financial circumstances, investment or trading objectives, risk tolerance, and liquidity needs. You are responsible for conducting your own independent research and seeking the advice of a qualified professional before making any investment or trading decisions.","category":"page"}]
}
