module VLQuantitativeFinancePackage

    # include -
    include("Include.jl")

    # export abstract types -
    export AbstractAssetModel
    export AbstractEquityPriceTreeModel
    export AbstractInterestRateTreeModel
    export AbstractContractModel
    export AbstractTreasuryDebtSecurity
    export AbstractCompoundingModel
    export AbstractStochasticChoiceProblem

    # export concrete types -
    export MyCRRLatticeNodeModel, MyGeometricBrownianMotionEquityModel, MyMultipleAssetGeometricBrownianMotionEquityModel, MyAdjacencyBasedCRREquityPriceTree, MyLongstaffSchwartzContractPricingModel, MyBlackScholesContractPricingModel
    export MyEuropeanCallContractModel, MyEuropeanPutContractModel, MyAmericanPutContractModel, MyAmericanCallContractModel, MyEquityModel
    export MyUSTreasuryZeroCouponBondModel, MyUSTreasuryCouponSecurityModel, DiscreteCompoundingModel, ContinuousCompoundingModel
    export MySymmetricBinaryInterestRateLatticeModel, MyBinaryInterestRateLatticeNodeModel
    export MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem, MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem
    export MyBiomialLatticeEquityNodeModel, MyBinomialEquityPriceTree
    export MySingleIndexModel, AbstractReturnModel
    export RealWorldBinomialProbabilityMeasure, RiskNeutralBinomialProbabilityMeasure, AbstractProbabilityMeasure

    # export functions/methods
    export build, payoff, profit, premium, sample, sample_endpoint, price, strip, populate, solve, YTM
    export estimate_implied_volatility
    
    # export the greeks -
    export theta, delta, gamma, vega, rho
end
