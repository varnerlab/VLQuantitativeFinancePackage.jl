# setup my paths -
const _PATH_TO_SRC = dirname(pathof(@__MODULE__));

# load external dependencies -
using Distributions
using Statistics
using LinearAlgebra
using DataFrames
using JLD2
using CSV
using MadNLP
using JuMP
using LsqFit
using Optim
using DataStructures
using StatsBase

# load my codes -
include(joinpath(_PATH_TO_SRC, "Types.jl"));
include(joinpath(_PATH_TO_SRC, "Factory.jl"));
include(joinpath(_PATH_TO_SRC, "Compute.jl"));
include(joinpath(_PATH_TO_SRC, "Longstaff.jl"));
include(joinpath(_PATH_TO_SRC, "Solve.jl"));
include(joinpath(_PATH_TO_SRC, "Greeks.jl"));
include(joinpath(_PATH_TO_SRC, "YTM.jl"));
include(joinpath(_PATH_TO_SRC, "Volatility.jl"));
include(joinpath(_PATH_TO_SRC, "EulerMaruyamaMethod.jl"));
include(joinpath(_PATH_TO_SRC, "Hippo.jl"));
include(joinpath(_PATH_TO_SRC, "Base.jl"));
include(joinpath(_PATH_TO_SRC, "Bandits.jl"));
include(joinpath(_PATH_TO_SRC, "Wolfram.jl"));
include(joinpath(_PATH_TO_SRC, "RL.jl"));