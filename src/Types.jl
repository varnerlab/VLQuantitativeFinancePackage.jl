abstract type AbstractAssetModel end
abstract type AbstractEquityPriceTreeModel <: AbstractAssetModel end
abstract type AbstractInterestRateTreeModel <: AbstractAssetModel end
abstract type AbstractContractModel <: AbstractAssetModel end

# --- Equity models ------------------------------------------------------------------------ #
struct MyLocalExpectationRegressionModel 
    
    # parameters -
    a0::Float64
    a1::Float64
    a2::Float64
    a3::Float64
    a4::Float64

    function MyLocalExpectationRegressionModel(a0,a1,a2,a3,a4)
        this = new(a0,a1,a2,a3,a4)
    end
end

mutable struct MyCRRLatticeNodeModel

    # data -
    price::Float64
    probability::Float64
    intrinsic::Union{Nothing,Float64}
    extrinsic::Union{Nothing,Float64}

    # constructor -
    MyCRRLatticeNodeModel() = new();
end

mutable struct MyGeometricBrownianMotionEquityModel <: AbstractAssetModel

    # data -
    μ::Float64
    σ::Float64

    # constructor -
    MyGeometricBrownianMotionEquityModel() = new()
end

mutable struct MyAdjacencyBasedCRREquityPriceTree <: AbstractEquityPriceTreeModel

    # data -
    data::Union{Nothing, Dict{Int, MyCRRLatticeNodeModel}}
    connectivity::Union{Nothing, Dict{Int64, Array{Int64,1}}}
    levels::Union{Nothing, Dict{Int64,Array{Int64,1}}}
    u::Float64
    d::Float64
    p::Float64
    ΔT::Float64
    μ::Float64

    # constructor 
    MyAdjacencyBasedCRREquityPriceTree() = new()
end

mutable struct MyLongstaffSchwartzContractPricingModel <: AbstractEquityPriceTreeModel

    # data -
    S::Array{Float64,2}
    r̄::Float64
    Δt::Float64

    # constructor -
    MyLongstaffSchwartzContractPricingModel() = new();
end

mutable struct MyBlackScholesContractPricingModel <: AbstractAssetModel

    # data -
    r::Float64
    Sₒ::Float64

    # constructor -
    MyBlackScholesContractPricingModel() = new();
end
# ------------------------------------------------------------------------------------------- #

# --- Contract models ----------------------------------------------------------------------- #
mutable struct MyEuropeanCallContractModel <: AbstractContractModel

    # data -
    K::Float64
    sense::Union{Nothing, Int64}
    DTE::Union{Nothing,Float64}
    IV::Union{Nothing, Float64}
    premium::Union{Nothing, Float64}
    ticker::Union{Nothing,String}
    copy::Union{Nothing, Int64}

    # constructor -
    MyEuropeanCallContractModel() = new()
end

mutable struct MyEuropeanPutContractModel <: AbstractContractModel

    # data -
    K::Float64
    sense::Union{Nothing, Int64}
    DTE::Union{Nothing,Float64}
    IV::Union{Nothing, Float64}
    premium::Union{Nothing, Float64}
    ticker::Union{Nothing,String}
    copy::Union{Nothing, Int64}

    # constructor -
    MyEuropeanPutContractModel() = new()
end

mutable struct MyAmericanCallContractModel <: AbstractContractModel

    # data -
    K::Float64
    sense::Union{Nothing, Int64}
    DTE::Union{Nothing,Float64}
    IV::Union{Nothing, Float64}
    premium::Union{Nothing, Float64}
    ticker::Union{Nothing,String}
    copy::Union{Nothing, Int64}

    # constructor -
    MyAmericanCallContractModel() = new()
end

mutable struct MyAmericanPutContractModel <: AbstractContractModel

    # data -
    K::Float64
    sense::Union{Nothing, Int64}
    DTE::Union{Nothing,Float64}
    IV::Union{Nothing, Float64}
    premium::Union{Nothing, Float64}
    ticker::Union{Nothing,String}
    copy::Union{Nothing, Int64}

    # constructor -
    MyAmericanPutContractModel() = new()
end

mutable struct MyEquityModel <: AbstractAssetModel

    # data -
    ticker::String
    purchase_price::Float64
    current_price::Float64
    direction::Int64
    number_of_shares::Int64

    # constructor -
    MyEquityModel() = new()
end
# -------------------------------------------------------------------------------------------- #

# --- Term structure of interest rates  ------------------------------------------------------ #
# abstract types -
abstract type AbstractLatticeModel end
abstract type AbstractLatticeNodeModel end


# concrete types -

# nodes
mutable struct MyBinaryInterestRateLatticeNodeModel

    # data -
    probability::Float64
    rate::Float64
    price::Float64

    # constructor -
    MyBinaryInterestRateLatticeNodeModel() = new();
end

# tree
mutable struct MySymmetricBinaryInterestRateLatticeModel <: AbstractInterestRateTreeModel
    
    # data -
    u::Float64      # up-factor
    d::Float64      # down-factor
    p::Float64      # probability of an up move
    rₒ::Float64     # root value
    T::Int64        # number of levels in the tree (zero based)
    
    connectivity::Union{Nothing,Dict{Int64, Array{Int64,1}}}    # holds the connectivity of the tree
    levels::Union{Nothing,Dict{Int64,Array{Int64,1}}} # nodes on each level of the tree
    data::Union{Nothing, Dict{Int64, MyBinaryInterestRateLatticeNodeModel}} # holds data in the tree

    # constructor -
    MySymmetricBinaryInterestRateLatticeModel() = new()
end
# -------------------------------------------------------------------------------------------- #