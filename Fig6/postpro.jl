using JLD
using Plots
using LaTeXStrings
using LinearAlgebra: det
using CovarianceEstimation
using StatsBase: countmap

for res_number in 1:21

    ################################################################################

    N_range         = 1:100
    p_range         = LinRange(0.05,0.95,50)
    q_range         = LinRange(0.0,0.2,50)
    sigma_range     = LinRange(0.0,0.1,50)
    tauD_range      = LinRange(0.0,0.2,50)

    N_values = Float64[]
    p_values = Float64[]
    q_values = Float64[]
    sigma_values = Float64[]
    tau_values = Float64[]

    L = -Inf
    N_MAP = 0
    p_MAP = 0.0
    q_MAP = 0.0
    sigma_MAP = 0.0
    tau_MAP = 0.0

    for i in 1:30
        obj = load(string("res1\\",i,"_",res_number,".jld"))
        res = obj["res"]
        N_values = append!(N_values, res[1])
        p_values = append!(p_values, res[2])
        q_values = append!(q_values, res[3])
        sigma_values = append!(sigma_values, res[4])
        tau_values = append!(tau_values, res[5])
        if res[11]>L
            N_MAP = res[6]
            p_MAP = res[7]
            q_MAP = res[8]
            sigma_MAP = res[9]
            tau_MAP = res[10]
            L = res[11]
        end
    end



    N_posterior    = zeros(length(N_range))
    p_posterior     = zeros(length(p_range))
    q_posterior     = zeros(length(q_range))
    sigma_posterior = zeros(length(sigma_range))
    tauD_posterior  = zeros(length(tauD_range))

    N_res = countmap(N_values)
    p_res = countmap(p_values)
    q_res = countmap(q_values)
    sigma_res = countmap(sigma_values)
    tauD_res = countmap(tau_values)

    for i in 1:length(N_range)
        N_posterior[i] = get(N_res,N_range[i],0)
    end
    for i in 1:length(p_range)
        p_posterior[i] = get(p_res,p_range[i],0)
    end
    for i in 1:length(q_range)
        q_posterior[i] = get(q_res,q_range[i],0)
    end
    for i in 1:length(sigma_range)
        sigma_posterior[i] = get(sigma_res,sigma_range[i],0)
    end
    for i in 1:length(tauD_range)
        tauD_posterior[i] = get(tauD_res,tauD_range[i],0)
    end

    ################################################################################
    default(palette = palette(:default))
    p2 = plot(N_range, N_posterior/sum(N_posterior), yaxis=nothing,ylabel=L"p(N)",legend=false,linewidth = 2,
    xlabel=L"N[-]")
    p3 = plot(p_range, p_posterior/sum(p_posterior),ylabel=L"p(p)", yaxis=nothing,legend=false,xticks=[0.2,0.5,0.8],linewidth = 2,
    xlabel=L"p[-]")
    p4 = plot(q_range, q_posterior/sum(q_posterior),ylabel=L"p(q)",legend=false, yaxis=nothing,linewidth = 2,xticks=[0.0,1e-11],
    xlabel=L"q[A]")
    p5 = plot(sigma_range, sigma_posterior/sum(sigma_posterior),ylabel=L"p(\sigma)",legend=false, yaxis=nothing,linewidth = 2,xticks=[0.0,3e-12],
    xlabel=L"\sigma[A]")
    p6 = plot(tauD_range, tauD_posterior/sum(tauD_posterior),ylabel=L"p(\tau_D)", yaxis=nothing,xticks=[0.0,0.1,0.2],linewidth = 2,legend=false,
    xlabel=L"\tau_D[s]")
    display(plot(p2, p3, p4, p5, p6, layout = 6,dpi=300,size=(400,350)))

    method = LinearShrinkage(DiagonalUnequalVariance(), 0.5)
    samples = [N_values';p_values';q_values';sigma_values';tau_values']


    covar = cov(method,samples')
    ent = 0.5*log(det(2*pi*ℯ*covar))
    print(ent)
    print('\n')

end