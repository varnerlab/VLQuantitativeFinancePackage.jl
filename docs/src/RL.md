# Reinforcement Learning
Fill me in

## Grid World
```@docs
VLQuantitativeFinancePackage.MyPeriodicRectangularGridWorldModel
VLQuantitativeFinancePackage.build(model::Type{MyPeriodicRectangularGridWorldModel}, data::NamedTuple)
```

## Wolfram policies and grids
Fill me in

```@docs
VLQuantitativeFinancePackage.MyOneDimensionalElementarWolframRuleModel
VLQuantitativeFinancePackage.build(model::Type{MyOneDimensionalElementarWolframRuleModel}, data::NamedTuple)
VLQuantitativeFinancePackage.MyTwoDimensionalTotalisticWolframRuleModel
VLQuantitativeFinancePackage.build(model::Type{MyTwoDimensionalTotalisticWolframRuleModel}, data::NamedTuple)
VLQuantitativeFinancePackage.solve(rulemodel::MyTwoDimensionalTotalisticWolframRuleModel, initialstate::Array{Int64,2}; steps::Int64 = 100)
```

## Wolfram Q-learning
Fill me in

```@docs
VLQuantitativeFinancePackage.MyWolframRuleQLearningAgentModel
VLQuantitativeFinancePackage.build(model::Type{MyWolframRuleQLearningAgentModel}, data::NamedTuple)
VLQuantitativeFinancePackage.MyWolframGridWorldModel
VLQuantitativeFinancePackage.build(model::Type{MyWolframGridWorldModel}, data::NamedTuple)
VLQuantitativeFinancePackage.sample(agent::MyWolframRuleQLearningAgentModel, environment::MyWolframGridWorldModel; maxsteps::Int = 100,
    ϵ::Float64 = 0.2)
```