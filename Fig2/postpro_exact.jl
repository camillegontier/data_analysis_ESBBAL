using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 400 #Number of data points
M = 187 #Number of simulations

# Ground truth parameters

N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

# Myopic protocol #########################################################
tau_entropy_values_Myopic_exact = zeros(T,M)
entropy_values_Myopic_exact = zeros(T,M)

N_MAP_values_Myopic_exact = zeros(T,M)
p_MAP_values_Myopic_exact = zeros(T,M)
q_MAP_values_Myopic_exact = zeros(T,M)
sigma_MAP_values_Myopic_exact = zeros(T,M)
tau_MAP_values_Myopic_exact = zeros(T,M)

computation_time_values_Myopic_exact = zeros(T,M)
ISI_values_Myopic_exact = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\Myopic_exact_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][7])
    for i in 1:T
        tau_entropy_values_Myopic_exact[i,j] = data[i][1]

        N_MAP_values_Myopic_exact[i,j] = data[i][2]
        p_MAP_values_Myopic_exact[i,j] = data[i][3]
        q_MAP_values_Myopic_exact[i,j] = data[i][4]
        sigma_MAP_values_Myopic_exact[i,j] = data[i][5]
        tau_MAP_values_Myopic_exact[i,j] = data[i][6]

        ISI_values_Myopic_exact[i,j] = ISI[i]
        computation_time_values_Myopic_exact[i,j] = data[i][8]

        entropy_values_Myopic_exact[i,j] = data[i][9]
    end
end

tau_entropy_Myopic_exact = mean(tau_entropy_values_Myopic_exact,dims=2)
entropy_Myopic_exact = mean(entropy_values_Myopic_exact,dims=2)

N_MAP_Myopic_exact = mean(N_MAP_values_Myopic_exact,dims=2)
p_MAP_Myopic_exact = mean(p_MAP_values_Myopic_exact,dims=2)
q_MAP_Myopic_exact = mean(q_MAP_values_Myopic_exact,dims=2)
sigma_MAP_Myopic_exact = mean(sigma_MAP_values_Myopic_exact,dims=2)
tau_MAP_Myopic_exact = mean(tau_MAP_values_Myopic_exact,dims=2)

###############################################################################
# @save "data_exact.jld"
@load "data.jld"
@load "data_exact.jld"

################################################################################

p6 = plot(entropy_constant_joint[1:200],label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300, title="Joint distribution",
ribbon=std(entropy_values_constant_joint[1:200,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Uniform_joint[1:200],label="Best Uniform",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_Uniform_joint[1:200,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_exponential_joint[1:200],label="Best Exponential",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_exponential_joint[1:200,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic[1:200],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,
color="black",line = :dash,
ribbon=std(entropy_values_Myopic[1:200,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic_exact[1:200],label="ESB-BAL (\"exact\")",linewidth = 2, thickness_scaling = 1,
color="grey",line = :dashdot,
ribbon=std(entropy_values_Myopic_exact[1:200,:],dims=2)/sqrt(M),fillalpha=.5)
savefig("exact.png")
