"""
    function solve(model::MySharpeRatioPortfolioChoiceProblem)::Dict{String,Any}

This method solves the Sharpe ratio maximization portfolio choice problem for a given instance of the [`MySharpeRatioPortfolioChoiceProblem`](@ref) model type.
This problem is encoded as a second-order cone program (SOCP) and solved using the COSMO solver via the JuMP modeling interface.

### Arguments
- `model::MySharpeRatioPortfolioChoiceProblem`: An instance of the [`MySharpeRatioPortfolioChoiceProblem`](@ref) that defines the problem parameters.

### Returns
- `Dict{String, Any}`: A dictionary containing the optimization results.

The results dictionary has the following keys:
- `"sharpe_ratio"`: The maximum Sharpe ratio achieved by the optimal portfolio.
- `"argmax"`: The optimal portfolio weights that achieve the maximum Sharpe ratio.
- `"numerator"`: The numerator of the Sharpe ratio at the optimal solution (excess return over risk-free rate).
- `"denominator"`: The denominator of the Sharpe ratio at the optimal solution (standard deviation of returns).
- `"status"`: The status of the optimization
"""
function solve(model::MySharpeRatioPortfolioChoiceProblem)::Dict{String,Any}
    
    # initialize -
    results = Dict{String,Any}()
    Σ = model.Σ;
    rfr = model.risk_free_rate;
    α = model.alpha;
    β = model.beta;
    gₘ = model.gₘ;

    # setup the problem -
    d = length(α);
    c = α .+ β .* gₘ .- rfr .* ones(d);

    # Cholesky: Sigma = U' * U  (Julia's cholesky gives Upper U)
    U = cholesky(Symmetric(Σ)).U;

    # setup the optimization model -
    opt_model = Model(()->COSMO.Optimizer(verbose=0, max_iter=500))
    @variable(opt_model, w[1:d] >= 0)
    @constraint(opt_model, sum(w) == 1)
    
    # SOC: ||U*w||_2 ≤ 1
    @constraint(opt_model, [1.0; U * w] in SecondOrderCone())
    @objective(opt_model, Max, dot(c, w))
    optimize!(opt_model)

    # check: was the optimization successful?
    @assert is_solved_and_feasible(opt_model)

    # get some results -
    w_opt = value.(w)                  # long-only, budgeted weights
    sr_num = dot(c, w_opt)             # numerator
    sr_den = norm(U * w_opt)           # ≈ 1 at optimum (up to tolerances)
    sr_opt = sr_num / sr_den           # Sharpe at the solution

    # package results -
    results["sharpe_ratio"] = sr_opt;
    results["argmax"] = w_opt;
    results["numerator"] = sr_num;
    results["denominator"] = sr_den;
    results["status"] = termination_status(opt_model);

    # return -
    return results;
end