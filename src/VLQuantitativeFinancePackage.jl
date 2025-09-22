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
    export AbstractMarkovModel
    export AbstractSamplingModel
    export AbstractWorldModel
    export AbstractPolicyModel
    export AbstractLearningModel


    # export concrete types -
    export MyCRRLatticeNodeModel, MyGeometricBrownianMotionEquityModel, MyMultipleAssetGeometricBrownianMotionEquityModel, MyAdjacencyBasedCRREquityPriceTree, MyLongstaffSchwartzContractPricingModel, MyBlackScholesContractPricingModel
    export MyEuropeanCallContractModel, MyEuropeanPutContractModel, MyAmericanPutContractModel, MyAmericanCallContractModel, MyEquityModel
    export MyUSTreasuryZeroCouponBondModel, MyUSTreasuryCouponSecurityModel, DiscreteCompoundingModel, ContinuousCompoundingModel
    export MySymmetricBinaryInterestRateLatticeModel, MyBinaryInterestRateLatticeNodeModel
    export MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem, MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem
    export MyBiomialLatticeEquityNodeModel, MyBinomialEquityPriceTree
    export MySingleIndexModel, AbstractReturnModel
    export RealWorldBinomialProbabilityMeasure, RiskNeutralBinomialProbabilityMeasure, AbstractProbabilityMeasure
    export MyOrnsteinUhlenbeckModel, MyHestonModel, EulerMaruyamaMethod
    export MySisoLegSHippoModel, estimate_hippo_parameters, prediction
    export MyGeneralAdjacencyRecombiningCommodityPriceTree
    
    # Markov models, MDPs, Bandits types and methods -
    export MyHiddenMarkovModel,MyPeriodicRectangularGridWorldModel
    export MyEpsilonSamplingBanditModel, MyTickerPickerWorldModel, MyTickerPickerRiskAwareWorldModel
    export preference

    # wolfram rules etc
    export MyOneDimensionalTotalisticWolframRuleModel, MyOneDimensionalElementarWolframRuleModel
    export MyTwoDimensionalElementaryWolframRuleModel, MyTwoDimensionalTotalisticWolframRuleModel
    export MyWolframRuleQLearningAgentModel, MyWolframGridWorldModel
    
    # Base functions -
    export log_growth_matrix

    # export functions/methods
    export build, payoff, profit, premium, sample, sample_endpoint, price, strip, populate, solve, YTM, typicalprice, discount
    export estimate_implied_volatility
    export expectation, variance
    
    # export the greeks -
    export theta, delta, gamma, vega, rho

    # data functions -
    export MyTrainingMarketDataSet;
end
