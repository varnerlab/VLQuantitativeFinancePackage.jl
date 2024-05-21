function _hippo_objective_function(p, signal, hidden);

    # run the model -
    (number_of_time_steps, number_of_hidden_states) = size(hidden);
    Ŷ = zeros(number_of_time_steps);
    
    # compute the output -
    for i ∈ 2:number_of_time_steps
        Ŷ[i] = dot(p, hidden[i,:]);
    end

    # compute the error -
    error_term = sum((signal - Ŷ).^2);

    # return -
    return error_term;
end

function solve(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64})::Tuple

    # initialize -
    Â = model.Â
    B̂ = model.B̂
    Ĉ = model.Ĉ
    D̂ = model.D̂
    Xₒ = model.Xₒ
    number_of_hidden_states = model.n;

    # build the time array -
    tₒ = tspan.start;
    tₙ = tspan.stop;
    dt = tspan.step;
    T = range(tₒ, step=dt, stop=tₙ) |> collect

    # initialize the state and output arrays -
    number_of_time_steps = length(T);
    Y = zeros(number_of_time_steps);
    X = zeros(number_of_time_steps, number_of_hidden_states);

    for i ∈ 1:number_of_hidden_states
        X[1,i] = Xₒ[i];
    end

    # main loop -
    for i ∈ 2:number_of_time_steps
        X[i,:] = Â*X[i-1,:]+B̂*signal[i-1];
        Y[i] = dot(Ĉ, X[i,:]);
    end

    # return the time and state arrays -
    return (T, X, Y);
end

function estimate_hippo_parameters(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64};
    method = LBFGS())
    
    # initialize -
    p = model.Ĉ;
 
    # solve the model to get the initial guess -
    (T, X, Y) = solve(model, tspan, signal);

    # setup the objective function -
    loss(p) = _hippo_objective_function(p,signal,X);
 
    # call the optimizer -
    opt_result = Optim.optimize(loss, p, method);
 
    # return the result -
    return Optim.minimizer(opt_result);
end