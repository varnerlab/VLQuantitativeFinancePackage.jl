function _build(modeltype::Type{T}, data::NamedTuple) where T <: Union{AbstractEquityPriceTreeModel, AbstractAssetModel, AbstractTreasuryDebtSecurity, AbstractStochasticChoiceProblem, AbstractReturnModel}
    
    # build an empty model
    model = modeltype();

    # if we have options, add them to the contract model -
    if (isempty(data) == false)
        for key ∈ fieldnames(modeltype)
            
            # check the for the key - if we have it, then grab this value
            value = nothing
            if (haskey(data, key) == true)
                # get the value -
                value = data[key]
            end

            # set -
            setproperty!(model, key, value)
        end
    end
 
    # return -
    return model
end

"""
    _build_nodes_level_dictionary(levels::Int64) -> Dict{Int64,Array{Int64,1}}
"""
function _build_nodes_level_dictionary(levels::Int64)::Dict{Int64,Array{Int64,1}}

    # initialize -
    index_dict = Dict{Int64, Array{Int64,1}}()

    counter = 0
    for l = 0:levels
        
        # create index set for this level -
        index_array = Array{Int64,1}()
        for _ = 1:(l+1)
            counter = counter + 1
            push!(index_array, counter)
        end

        index_dict[l] = (index_array .- 1) # zero based
    end

    # return -
    return index_dict
end

function _build_connectivity_dictionary(h::Int)::Dict{Int64, Array{Int64,1}}

    # compute connectivity - 
    number_items_per_level = [i for i = 1:(h+1)]
    tmp_array = Array{Int64,1}()
    theta = 0
    for value in number_items_per_level
        for _ = 1:value
            push!(tmp_array, theta)
        end
        theta = theta + 1
    end

    N = sum(number_items_per_level[1:(h)])
    connectivity_index_array = Array{Int64,2}(undef, N, 3)
    for row_index = 1:N

        # index_array[row_index,1] = tmp_array[row_index]
        connectivity_index_array[row_index, 1] = row_index
        connectivity_index_array[row_index, 2] = row_index + 1 + tmp_array[row_index]
        connectivity_index_array[row_index, 3] = row_index + 2 + tmp_array[row_index]
    end
    
    # adjust for zero base -
    zero_based_array = connectivity_index_array .- 1;

    # build connectivity dictionary -
    N = sum(number_items_per_level[1:end-1])
    connectivity = Dict{Int64, Array{Int64,1}}()
    for i ∈ 0:(N-1)
        # grab the connectivity -
        connectivity[i] = reverse(zero_based_array[i+1,2:end])
    end

    # put it back in order -
    for i ∈ 0:(N-1)
        # grab the connectivity -
        connectivity[i] = zero_based_array[i+1,2:end]
    end

    return connectivity
end

# Old IMPL that works - do not delete, just in case .... switched to two phase impl, *seems* to be ok??
# However, we need this here for backward compatibility
function build(modeltype::Type{MyAdjacencyBasedCRREquityPriceTree}; 
    h::Int = 1, μ::Float64 = 0.01, σ::Float64 = 0.1, T::Float64 = (1.0/365.0), 
    Sₒ::Float64 = 1.0)::MyAdjacencyBasedCRREquityPriceTree
     
    # initialize -
    model = MyAdjacencyBasedCRREquityPriceTree(); # this model is empty
    nodes_dictionary = Dict{Int, MyCRRLatticeNodeModel}()

    # compute u, d and p
    ΔT = T / h
    u = exp(σ * sqrt(ΔT))
    d = 1.0/u;
    p = (exp(µ * ΔT) - d) / (u - d)

    @show (ΔT,u,d,p)
  
    # # compute connectivity - 
    # number_items_per_level = [i for i = 1:(h+1)]
    # tmp_array = Array{Int64,1}()
    # theta = 0
    # for value in number_items_per_level
    #     for _ = 1:value
    #         push!(tmp_array, theta)
    #     end
    #     theta = theta + 1
    # end

    # N = sum(number_items_per_level[1:(h)])
    # connectivity_index_array = Array{Int64,2}(undef, N, 3)
    # for row_index = 1:N

    #     # index_array[row_index,1] = tmp_array[row_index]
    #     connectivity_index_array[row_index, 1] = row_index
    #     connectivity_index_array[row_index, 2] = row_index + 1 + tmp_array[row_index]
    #     connectivity_index_array[row_index, 3] = row_index + 2 + tmp_array[row_index]
    # end
    
    # # adjust for zero base -
    # zero_based_array = connectivity_index_array .- 1;

    # # build connectivity dictionary -
    # N = sum(number_items_per_level[1:end-1])
    # connectivity = Dict{Int64, Array{Int64,1}}()
    # for i ∈ 0:(N-1)
    #     # grab the connectivity -
    #     connectivity[i] = reverse(zero_based_array[i+1,2:end])
    # end

    # compute the price and probability, and store in the nodes dictionary
    counter = 0;
    for t ∈ 0:h

        # prices -
        for k ∈ 0:t
            
            t′ = big(t)
            k′ = big(k)

            # compute the prices and P for this level
            price = Sₒ*(u^(t-k))*(d^(k));
            P = binomial(t′,k′)*(p^(t-k))*(1-p)^(k);

            # create a node model -
            node = MyCRRLatticeNodeModel();
            node.price = price
            node.probability = P;
            node.intrinsic = 0.0; # intrinsic value gets updated later, for now -> 0.0
            node.extrinsic = 0.0; # extrinsic value gets updated later, for now -> 0.0
            
            # push this into the array -
            nodes_dictionary[counter] = node;
            counter += 1
        end
    end

    # # put it back in order -
    # for i ∈ 0:(N-1)
    #     # grab the connectivity -
    #     connectivity[i] = zero_based_array[i+1,2:end]
    # end

    # # set the data, and connectivity for the model -
    model.data = nodes_dictionary;
    # model.connectivity = connectivity;
    model.connectivity = _build_connectivity_dictionary(h)
    model.levels = _build_nodes_level_dictionary(h)
    model.p = p;
    model.u = u;
    model.ΔT = ΔT
    model.μ = μ
    model.d = d
    model.T = T

    # return -
    return model
end

# short build methods -
build(model::Type{MyEuropeanCallContractModel}, data::NamedTuple)::MyEuropeanCallContractModel = _build(model, data)
build(model::Type{MyEuropeanPutContractModel}, data::NamedTuple)::MyEuropeanPutContractModel = _build(model, data)
build(model::Type{MyAmericanPutContractModel}, data::NamedTuple)::MyAmericanPutContractModel = _build(model, data)
build(model::Type{MyAmericanCallContractModel}, data::NamedTuple)::MyAmericanCallContractModel = _build(model, data)
build(model::Type{MyLongstaffSchwartzContractPricingModel}, data::NamedTuple)::MyLongstaffSchwartzContractPricingModel = _build(model, data)
build(model::Type{MyGeometricBrownianMotionEquityModel}, data::NamedTuple)::MyGeometricBrownianMotionEquityModel = _build(model, data)
build(model::Type{MyMultipleAssetGeometricBrownianMotionEquityModel}, data::NamedTuple)::MyMultipleAssetGeometricBrownianMotionEquityModel = _build(model, data)
build(model::Type{MyBlackScholesContractPricingModel}, data::NamedTuple)::MyBlackScholesContractPricingModel = _build(model, data)
build(model::Type{MySymmetricBinaryInterestRateLatticeModel}, data::NamedTuple)::MySymmetricBinaryInterestRateLatticeModel = _build(model, data);
build(model::Type{MyBinaryInterestRateLatticeNodeModel}, data::NamedTuple)::MyBinaryInterestRateLatticeNodeModel = _build(model, data);
build(model::Type{MyUSTreasuryCouponSecurityModel}, data::NamedTuple)::MyUSTreasuryCouponSecurityModel = _build(model, data);
build(model::Type{MyUSTreasuryZeroCouponBondModel}, data::NamedTuple)::MyUSTreasuryZeroCouponBondModel = _build(model, data);
build(model::Type{MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem}, data::NamedTuple)::MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem = _build(model, data);
build(model::Type{MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem}, data::NamedTuple)::MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem = _build(model, data);
build(model::Type{MyAdjacencyBasedCRREquityPriceTree}, data::NamedTuple)::MyAdjacencyBasedCRREquityPriceTree = _build(model, data);
build(model::Type{MyBinomialEquityPriceTree}, data::NamedTuple)::MyBinomialEquityPriceTree = _build(model, data);
build(model::Type{MySingleIndexModel}, dataL::NamedTuple) = _build(model, data);
build(model::Type{MyHestonModel}, data::NamedTuple)::MyHestonModel = _build(model, data);


function build(modeltype::Type{MyOrnsteinUhlenbeckModel}, data::NamedTuple)::MyOrnsteinUhlenbeckModel

    # initialize -
    model = modeltype();

    model.μ = data.μ
    model.σ = data.σ
    model.θ = data.θ

    # return -
    return model;
end

function build(modeltype::Type{MySisoLegSHippoModel}, data::NamedTuple)::MySisoLegSHippoModel

    # initialize -
    model = modeltype(); # build an empty model

    # get data -
    number_of_hidden_states = data.number_of_hidden_states;
    Δt = data.Δt;
    uₒ = data.uₒ;
    C = data.C;

    # A matrix -
    A = zeros(number_of_hidden_states,number_of_hidden_states);
    for i ∈ 1:number_of_hidden_states
        for j ∈ 1:number_of_hidden_states
        
            a = -sqrt((2*i+1))*sqrt((2*j+1));
            b = nothing;
            if (i > j)
                b = 1;
            elseif (i == j)
                b = (i+1)/(2*i+1);
            else
                b = 0
            end
            A[i,j] = a*b;
        end
    end

    # B matrix -
    B = zeros(number_of_hidden_states);
    for i ∈ 1:number_of_hidden_states
        B[i] = (2*i+1) |> sqrt
    end


    # discretize the arrays using the Bilinear method -
    Â = inv((I - (Δt/2)*A))*(I + (Δt/2)*A);
    B̂ = inv((I - (Δt/2)*A))*(Δt)*B;
    Ĉ = C; # initialize a random C matrix (user can update this later)
    D̂ = zeros(number_of_hidden_states); # initialize a zero D matrix (user can update this later)

    # set the values -
    model.Â = Â;
    model.B̂ = B̂;
    model.Ĉ = Ĉ;
    model.D̂ = D̂;
    model.n = number_of_hidden_states;
    model.Xₒ = B̂*uₒ;

    # return -
    return model;
end