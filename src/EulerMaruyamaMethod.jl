
function _solve(model::MyHestonModel, tspan::NamedTuple,
    initialconditions::AbstractArray, N::Int64, method::EulerMaruyamaMethod)::Tuple
    
    return (0,0);
end


function _solve(model::MyOrnsteinUhlenbeckModel, tspan::NamedTuple,
    initialconditions::AbstractArray, N::Int64, method::EulerMaruyamaMethod)::Tuple

    # initialize -
    μ = model.μ
    σ = model.σ
    θ = model.θ
    Xₒ = initialconditions;

    # build the time array -
    tₒ = tspan.start;
    tₙ = tspan.stop;
    dt = tspan.step;
    T = range(tₒ, step=dt, stop=tₙ) |> collect
    M = length(T) # how many time steps do we have?
    X = Array{Float64,2}(undef, M, N) # initialize an empty array to store the price paths

    # fill in the first row, this is the initial price
    foreach(p -> X[1,p] = Xₒ[1], 1:N)

    # pre generate the noise matrix -
    ZM = Normal(0,1) |> d->rand(d, M-1, N)
    
    # update the state array -
    for p ∈ 1:N
        for i ∈ 2:M
            X[i,p] = X[i-1,p] + θ*(μ - X[i-1,p])*dt + σ*sqrt(dt)*ZM[i,p] # update
        end
    end

    # return the time and state arrays -
    return (T, X);
end

function solve(model::AbstractAssetModel, tspan::NamedTuple,
    initialconditions::AbstractArray; method::AbstractStochasticSolverModel = EulerMaruyamaMethod(), 
    N::Int64 = 100)::Tuple
    return _solve(model, tspan, initialconditions, N, method);
end