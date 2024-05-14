function _hippo_objective_function(p, model, tspan, signal)

    # run the model -
    model.Ĉ = p;
    Ŷ = solve(model, tspan, signal)[3]; # get the predicted signal -

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
 
    # setup the objective function -
    loss(p) = _hippo_objective_function(p, model, tspan, signal);
 
    # call the optimizer -
    opt_result = Optim.optimize(loss, p, GradientDescent());
 
    # return the result -
    return Optim.minimizer(opt_result);
end