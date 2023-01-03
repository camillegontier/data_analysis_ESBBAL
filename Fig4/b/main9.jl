using Distributions
using StatsBase
using LaTeXStrings
using JLD
using DelimitedFiles

include("metropolis_hastings.jl")

T = 234

file_name = "cell3_det_100Hz_long.txt"
data = readdlm(joinpath("Data",file_name), ',', '\n')
delta_t = data[1:T,2]
e = -data[1:T,1]
e = e/maximum(e)

# Initialization #############################################################

# number of iterations of the algorithm.
nb_it   = 10000

# Discretized range of values for the posterior
N_range         = 1:100
p_range         = LinRange(0.05,0.95,50)
q_range         = LinRange(0.0,0.2,50)
sigma_range     = LinRange(0.0,0.1,50)
tauD_range      = LinRange(0.0,0.2,50)

# Runs MH #####################################################################

results = MH(e,N_range,p_range,q_range,sigma_range,tauD_range,delta_t,nb_it)

# Postpro #####################################################################

save(string(Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]),"_9.jld"), "res", results)
