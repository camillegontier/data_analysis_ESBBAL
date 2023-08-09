using BinomialSynapses
using JLD
using StatsBase
using Distributions: Exponential
using LinearAlgebra: det

N = 47
p = 0.27
q = 1.0
sigma = 0.5
tau = 0.17

candidates = LinRange(0.005,2,64)

sim = NestedFilterSimulation(
        N, p, q, sigma, tau,
        1:2:100,
        LinRange(0.05,0.95,50),
        LinRange(0.01,2,50),
        LinRange(0.01,1,50),
        LinRange(0.01,0.2,50),
        1024*4, 256*4, 40,
        timestep = Myopic(candidates,0.0,"MAP")
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

        map = MAP(sim.fstate.model)
        map_N = map.N
        map_p = map.p
        map_q = map.q
        map_σ = map.σ
        map_τ = map.τ

        s = [Nrng[Nind]'./Nrng[end];prng[pind]'./prng[end];qrng[qind]'./qrng[end];σrng[σind]'./σrng[end];τrng[τind]'./τrng[end]]
        Σ_est = cov(s')
        entropy_gaussian = 0.5*log(det(2*pi*ℯ*Σ_est))

        return map_N, map_p, map_q, map_σ, map_τ, sim.times, time, entropy_gaussian
end

function f2(data)
        save(string("Myopic_",Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"]),".jld"), "data", data)
end

rec = Recording(f1, f2, sim)

sim.hstate.n .= N
times, epsps = run!(sim, T = 200, recording = Recording(f1, f2, sim))