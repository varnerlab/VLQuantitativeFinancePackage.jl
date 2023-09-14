"""
    solve(problem::MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem) -> Dict{String,Any}
"""
function solve(problem::MyMarkowitzRiskyAssetOnlyPortfiolioChoiceProblem)::Dict{String,Any}

    # initialize -
    results = Dict{String,Any}()
    Σ = problem.Σ;
    μ = problem.μ;
    R = problem.R;
    bounds = problem.bounds;
    wₒ = problem.initial

    # setup the problem -
    d = length(μ)
    model = Model(()->MadNLP.Optimizer(print_level=MadNLP.ERROR, max_iter=500))
    @variable(model, bounds[i,1] <= w[i=1:d] <= bounds[i,2], start=wₒ[i])

    # set objective function -
    @objective(model, Min, transpose(w)*Σ*w);

    # setup the constraints -
    @constraints(model, 
        begin
            # my turn constraint
            transpose(μ)*w >= R
            sum(w) == 1.0
        end
    );

    # run the optimization -
    optimize!(model)

    # populate -
    w_opt = value.(w);
    results["argmax"] = w_opt
    results["reward"] = transpose(μ)*w_opt; 
    results["objective_value"] = objective_value(model);
    results["status"] = termination_status(model);

    # return -
    return results
end

"""
    solve(problem::MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem) -> Dict{String,Any}
"""
function solve(problem::MyMarkowitzRiskyRiskFreePortfiolioChoiceProblem)::Dict{String,Any}

    # initialize -
    results = Dict{String,Any}()
    Σ = problem.Σ;
    μ = problem.μ;
    R = problem.R;
    bounds = problem.bounds;
    initial = problem.initial
    rfr = problem.risk_free_rate

    # setup the problem -
    d = length(μ)
    model = Model(()->MadNLP.Optimizer(print_level=MadNLP.ERROR, max_iter=500))
    @variable(model, bounds[i,1] <= w[i=1:d] <= bounds[i,2], start=initial[i])

    # set objective function -
    @objective(model, Min, transpose(w)*Σ*w);

    # setup the constraints -
    @constraints(model, 
        begin
            # my turn constraint
            transpose(μ)*w + (1.0 - sum(w))*rfr >= R
        end
    );

    # run the optimization -
    optimize!(model)

    # populate -
    w_opt = value.(w);
    results["reward"] = transpose(μ)*w_opt + (1.0 - sum(w_opt))*rfr;
    results["argmax"] = w_opt;
    results["objective_value"] = objective_value(model);
    results["status"] = termination_status(model);

    # return -
    return results
end