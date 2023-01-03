using BinomialSynapses
using JLD
using StatsBase
using Distributions: Exponential
using LinearAlgebra: det

N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

sim = NestedFilterSimulation(
        N, p, q, sigma, tau,
        1:30,
        LinRange(0.05,0.95,30),
        LinRange(0.1,2,30),
        LinRange(0.05,1,30),
        LinRange(0.05,1,30),
        1024, 256, 12,
        timestep = FixedTimestep(1.0)
        )

function f1(sim, time)
    Nind = Array(sim.fstate.model.Nind)
    Nrng = Array(sim.fstate.model.Nrng)


    pind = Array(sim.fstate.model.pind)
    prng = Array(sim.fstate.model.prng)


    qind = Array(sim.fstate.model.qind)
    qrng = Array(sim.fstate.model.qrng)


    σind = Array(sim.fstate.model.σind)
    σrng = Array(sim.fstate.model.σrng)


    τind = Array(sim.fstate.model.τind)
    τrng = Array(sim.fstate.model.τrng)
    τ_posterior = zeros(length(τrng))
    for j in 1:length(τrng)
        τ_posterior[j] = count(i->(i==j),τind)
    end
    entropy_τ = entropy(τ_posterior/sum(τ_posterior))

    map = MAP(sim.fstate.model)
    map_N = map.N
    map_p = map.p
    map_q = map.q
    map_σ = map.σ
    map_τ = map.τ

    s = [Nrng[Nind]';prng[pind]';qrng[qind]';σrng[σind]';τrng[τind]']
    Σ_est = cov(s')
    entropy_gaussian = 0.5*log(det(2*pi*ℯ*Σ_est))

    return entropy_τ, map_N, map_p, map_q, map_σ, map_τ, sim.times, time, entropy_gaussian
end

function f2(data)
        save(string("constant_joint_",Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]),".jld"), "data", data)
end

rec = Recording(f1, f2, sim)

sim.hstate.n .= N
times, epsps = run!(sim, T = 400, recording = Recording(f1, f2, sim))
