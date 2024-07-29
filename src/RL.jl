
# -- PRIVATE METHODS BELOW HERE ------------------------------------------------------------------------------------------- #
function _world(model::MyWolframGridWorldModel, s::Int, a::Int)::Tuple{Int,Float64}

    # initialize -
    s′ = nothing
    r = nothing
    
    # get data from the model -
    coordinates = model.coordinates;
    moves = model.moves
    states = model.states;
    rewards = model.rewards;

    # where are we now?
    current_position = coordinates[s];

    # get the perturbation -
    Δ = moves[a];
    new_position = current_position .+ Δ

    # before we go on, have we "driven off the grid"?
    if (haskey(states, new_position) == true)

        # lookup the new state -
        s′ = states[new_position];
        r = rewards[s′];
    else
       
        # ok: so we are all the grid. Bounce us back to to the current_position, and charge a huge penalty 
        s′ = states[current_position];
        r = -1000000000000.0
    end

    # return -
    return (s′,r);
end

function _update(model::MyWolframRuleQLearningAgentModel, data::NamedTuple)::MyWolframRuleQLearningAgentModel

    # grab the s,a,reward and next state from the data tuple
    s = data[:s];
    a = data[:a];
    r = data[:r];
    s′ = data[:s′];
    
    # grab parameters from the model -
    γ, Q, α = model.γ, model.Q, model.α

    # use the update rule to update Q -
    Q[s,a] += α*(r+γ*maximum(Q[s′,:]) - Q[s,a])

    # return -
    return model;
end
# -- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------------------- #

# -- PUBLIC METHODS BELOW HERE -------------------------------------------------------------------------------------------- #
function sample(agent::MyWolframRuleQLearningAgentModel, environment::AbstractWorldModel, startstate::Tuple{Int,Int}, maxsteps::Int;
    ϵ::Float64 = 0.2)::MyWolframRuleQLearningAgentModel

    # initialize -
    s = environment.states[startstate]
    actions = agent.actions;
    number_of_actions = length(actions);

    # simulation loop -
    for _ ∈ 1:maxsteps

        a = nothing;
        if rand() < ϵ

            # we generate a random action
            a = rand(1:number_of_actions);
        else

            # ok: so we are in some state s, let's use our memory to suggest a new action
            Q = agent.Q;
            a = argmax(Q[s,:]);
        end

        # check the action -
        s′, r = nothing, nothing;
        current_position = environment.coordinates[s];
        new_position = current_position .+ environment.moves[a]
        if (haskey(environment.states, new_position) == true)

            # ask the world, what is my next state and reward from this (s,a)
            (s′,r) = environment(s,a)
        else
            s′ = s;
            r = -1000000000000.0;
        end
        
        # update my model -
        agent = agent((
            s = s, a = a, r = r, s′ = s′
        ));

        # move -
        s = s′;
    end

    # return -
    return agent
end

(model::MyWolframRuleQLearningAgentModel)(data::NamedTuple) = _update(model, data);
(model::MyWolframGridWorldModel)(s::Int,a::Int) = _world(model, s, a);
# -- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------------------- #
