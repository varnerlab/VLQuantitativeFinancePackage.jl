
# -- PRIVATE METHODS BELOW HERE ------------------------------------------------------------------------------------------- #
function _world(model::MyWolframGridWorldModel, t::Int, s::Int, a::Int)::Tuple{Int, Float64}

    # initialize -
    s′ = nothing
    r = nothing
    
    # grab the parameters from the model -
    dataset = model.data;

    # what is the state, action and reward?
    data = dataset[t];
    a′ = data[end]+1; # the last element is the next state (correct for zero)
    if (a′ == a)
        r = 1.0;
    else
        r = -1.0;
    end

    # jump to the next state -
    s′ = rand(1:model.number_of_states);

    # return -
    return (s′,r);
end

function _update(model::MyWolframRuleQLearningAgentModel, data::NamedTuple)::MyWolframRuleQLearningAgentModel

    # grab the s,a,reward and next state from the data tuple
    s = data[:s];
    a = data[:a];
    r = data[:r];
    s′ = data[:s′];

    @show s,a,r,s′
    
    # grab parameters from the model -
    γ, Q, α = model.γ, model.Q, model.α

    # use the update rule to update Q -
    Q[s,a] += α*(r+γ*maximum(Q[s′,:]) - Q[s,a])

    # return -
    return model;
end
# -- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------------------- #

# -- PUBLIC METHODS BELOW HERE -------------------------------------------------------------------------------------------- #
function sample(agent::MyWolframRuleQLearningAgentModel, environment::AbstractWorldModel; maxsteps::Int = 100,
    ϵ::Float64 = 0.2)::MyWolframRuleQLearningAgentModel

    # initialize -
    s = 1; # for now, we start at the first state
    actions = agent.actions;
    number_of_actions = length(actions);

    # simulation loop -
    for t ∈ 1:maxsteps

        a = nothing;
        if rand() < ϵ

            # we generate a random action
            a = rand(1:number_of_actions);
        else

            # ok: so we are in some state s, let's use our memory to suggest a new action
            Q = agent.Q;
            a = argmax(Q[s,:]);
        end

        # initialize -
        s′, r = nothing, nothing; 
        
        # ask the world, what is my next state and reward from this (s,a)
        (s′,r) = environment(t, s, a)
        
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
(model::MyWolframGridWorldModel)(t::Int, s::Int, a::Int) = _world(model, t, s, a);
# -- PRIVATE METHODS ABOVE HERE ------------------------------------------------------------------------------------------- #
