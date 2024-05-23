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


"""
    prediction(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64,1}; L::Int64 = 10) -> Tuple

The `prediction` function predicts the output of the single input single output (SISO) HiPPO model given a signal vector using the Leg-S parameterization and bilinear discretization. 
The function returns the time array, the hidden state array and the output array. 

### Arguments
- `model::MySisoLegSHippoModel`: A model struct that defines the HiPPO model, see [MySisoLegSHippoModel](@ref) for details on the model struct.
- `tspan::NamedTuple`: A named tuple that defines the time span for the simulation. The named tuple should have the fields `start`, `stop` and `step`.
- `signal::Array{Float64}`: An array of input signals to the model.
- `S::Int64`: Circular buffer size used by the prediction function. The default value is 10 steps. After `L` steps, the input signal is reset to the first signal value.

### Returns
- `Tuple`: A tuple of the time array `T`, hidden state array `X` and the output array `Y`.
"""
function prediction(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64,1};
    S::Int64 = 10, B::Float64 = 40.0, α::Float64 = 0.25, β::Float64 = 0.10)::Tuple
    
    # initialize -
    Â = model.Â
    B̂ = model.B̂
    Ĉ = model.Ĉ
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

    # update the initial state -
    for i ∈ 1:number_of_hidden_states
        X[1,i] = Xₒ[i];
    end

    # main loop -
    for i ∈ 2:number_of_time_steps
        
        # what index in the signal array should we use?
        j = ((i-2) % S) + 1;
        u = signal[j]; # get the input 
        # u = Y[i-1]; # get the input

        # update the state and output -
        X[i,:] = Â*X[i-1,:]+B̂*u;
        Y[i] = dot(Ĉ, X[i,:]);

        # ok, so we some stability issues here, let's try to fix it -
        if (abs(Y[i]) ≥ B*(1+α*randn()))
            
            # reset the hidden states -
            for k ∈ 1:number_of_hidden_states
                X[i,k] = Xₒ[k]*(1+β*randn()); # jump back to the last stable state
            end
        end
    end

    # return the time and state arrays -
    return (T, X, Y);
end

"""
    solve(model::MySisoLegSHippoModel, tspan::NamedTuple, signal::Array{Float64}) -> Tuple

The `solve` function solves the single input single output (SISO) HiPPO model using the Leg-S parameterization and 
bilinear discretization. The function returns the time array, the hidden state array and the output array.

### Arguments
- `model::MySisoLegSHippoModel`: A model struct that defines the HiPPO model, see [MySisoLegSHippoModel](@ref) for details on the model struct.
- `tspan::NamedTuple`: A named tuple that defines the time span for the simulation. The named tuple should have the fields `start`, `stop` and `step`.
- `signal::Array{Float64}`: An array of input signals to the model.

### Returns
- `Tuple`: A tuple of the time array `T`, hidden state array `X` and the output array `Y`.
"""
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