using BinomialSynapses
using JLD
using StatsBase
using Distributions: Exponential
using LinearAlgebra: det

N = 10
p = 0.85
q = 1.0
sigma = 0.2
tau = 0.2

candidates = LinRange(0.005,2,64)

sim = NestedFilterSimulation(
        N, p, q, sigma, tau,
        1:20,
        LinRange(0.05,0.95,25),
        LinRange(0.1,2,25),
        LinRange(0.05,1,25),
        LinRange(0.05,1,25),
        1024, 256, 12,
        timestep = Myopic_tau(candidates,0.2)
        )

function f1(sim, time)
        τind = Array(sim.fstate.model.τind)
        τrng = Array(sim.fstate.model.τrng)
        τ_posterior = zeros(length(τrng))
        for j in 1:length(τrng)
            τ_posterior[j] = count(i->(i==j),τind)
        end
        entropy_τ = entropy(τ_posterior/sum(τ_posterior))

        return entropy_τ, sim.times
end

function f2(data)
        save(string("main_02_",Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]),".jld"), "data", data)
end

rec = Recording(f1, f2, sim)

sim.hstate.n .= N
times, epsps = run!(sim, T = 200, recording = Recording(f1, f2, sim))
