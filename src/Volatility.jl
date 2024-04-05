
function _iv_objective_function(p, contract::T, Sₒ::Float64, h::Int64, r̄::Float64)::Float64 where {T<:AbstractContractModel}
    
    # set the volatility -
    IVᵢ = p[1]
    DTE = contract.DTE; # this is measured in years
    premium = contract.premium;

    # build a lattice model with the current parameters, populate it and compute the premium -
    treemodel = build(MyAdjacencyBasedCRREquityPriceTree, 
        (μ = r̄, T = DTE, σ = IVᵢ)) |> (x-> populate(x, Sₒ = Sₒ, h = h));
    
    # compute the premium -
    𝒫 = premium(contract, treemodel);

    # compute the error -
    error_term = (𝒫 - premium)^2;
    
    # return -
    return error_term;
end


function IV(contract::T; Sₒ::Float64 = 100.0, h::Int64 = 1, r̄::Float64 = 0.05)::Float64 where {T<:AbstractContractModel}
    
    # initialize -
    IV₀ = contract.IV;
    p = [IV₀]; # initial guess for the implied volatility

    # setup the objective function -
    OF(p) = _iv_objective_function(p, contract, Sₒ, h, r̄);

    # call the optimizer -
    opt_result = Optim.optimize(OF, IVₒ, BFGS())

    # return the result -
    return Optim.minimizer(opt_result)[1]
end