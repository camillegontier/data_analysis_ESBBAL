using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 400 #Number of data points
M = 1500 #Number of simulations

# Constant protocol #########################################################

tau_entropy_values_constant_joint = zeros(T,M)
entropy_values_constant_joint = zeros(T,M)

N_MAP_values_constant_joint = zeros(T,M)
p_MAP_values_constant_joint = zeros(T,M)
q_MAP_values_constant_joint = zeros(T,M)
sigma_MAP_values_constant_joint = zeros(T,M)
tau_MAP_values_constant_joint = zeros(T,M)

computation_time_values_constant_joint = zeros(T,M)
ISI_values_constant_joint = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\constant_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][7])
    for i in 1:T
        tau_entropy_values_constant_joint[i,j] = data[i][1]

        N_MAP_values_constant_joint[i,j] = data[i][2]
        p_MAP_values_constant_joint[i,j] = data[i][3]
        q_MAP_values_constant_joint[i,j] = data[i][4]
        sigma_MAP_values_constant_joint[i,j] = data[i][5]
        tau_MAP_values_constant_joint[i,j] = data[i][6]

        ISI_values_constant_joint[i,j] = ISI[i]
        computation_time_values_constant_joint[i,j] = data[i][8]

        entropy_values_constant_joint[i,j] = data[i][9]
    end
end

tau_entropy_constant_joint = mean(tau_entropy_values_constant_joint,dims=2)
entropy_constant_joint = mean(entropy_values_constant_joint,dims=2)

N_MAP_constant_joint = mean(N_MAP_values_constant_joint,dims=2)
p_MAP_constant_joint = mean(p_MAP_values_constant_joint,dims=2)
q_MAP_constant_joint = mean(q_MAP_values_constant_joint,dims=2)
sigma_MAP_constant_joint = mean(sigma_MAP_values_constant_joint,dims=2)
tau_MAP_constant_joint = mean(tau_MAP_values_constant_joint,dims=2)

# Exponential protocol #########################################################

tau_entropy_values_exponential_joint = zeros(T,M)
entropy_values_exponential_joint = zeros(T,M)

N_MAP_values_exponential_joint = zeros(T,M)
p_MAP_values_exponential_joint = zeros(T,M)
q_MAP_values_exponential_joint = zeros(T,M)
sigma_MAP_values_exponential_joint = zeros(T,M)
tau_MAP_values_exponential_joint = zeros(T,M)

computation_time_values_exponential_joint = zeros(T,M)
ISI_values_exponential_joint = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\exponential_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][7])
    for i in 1:T
        tau_entropy_values_exponential_joint[i,j] = data[i][1]

        N_MAP_values_exponential_joint[i,j] = data[i][2]
        p_MAP_values_exponential_joint[i,j] = data[i][3]
        q_MAP_values_exponential_joint[i,j] = data[i][4]
        sigma_MAP_values_exponential_joint[i,j] = data[i][5]
        tau_MAP_values_exponential_joint[i,j] = data[i][6]

        ISI_values_exponential_joint[i,j] = ISI[i]
        computation_time_values_exponential_joint[i,j] = data[i][8]

        entropy_values_exponential_joint[i,j] = data[i][9]
    end
end

tau_entropy_exponential_joint = mean(tau_entropy_values_exponential_joint,dims=2)
entropy_exponential_joint = mean(entropy_values_exponential_joint,dims=2)

N_MAP_exponential_joint = mean(N_MAP_values_exponential_joint,dims=2)
p_MAP_exponential_joint = mean(p_MAP_values_exponential_joint,dims=2)
q_MAP_exponential_joint = mean(q_MAP_values_exponential_joint,dims=2)
sigma_MAP_exponential_joint = mean(sigma_MAP_values_exponential_joint,dims=2)
tau_MAP_exponential_joint = mean(tau_MAP_values_exponential_joint,dims=2)

# Uniform protocol #########################################################
tau_entropy_values_Uniform_joint = zeros(T,M)
entropy_values_Uniform_joint = zeros(T,M)

N_MAP_values_Uniform_joint = zeros(T,M)
p_MAP_values_Uniform_joint = zeros(T,M)
q_MAP_values_Uniform_joint = zeros(T,M)
sigma_MAP_values_Uniform_joint = zeros(T,M)
tau_MAP_values_Uniform_joint = zeros(T,M)

computation_time_values_Uniform_joint = zeros(T,M)
ISI_values_Uniform_joint = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\Uniform_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][7])
    for i in 1:T
        tau_entropy_values_Uniform_joint[i,j] = data[i][1]

        N_MAP_values_Uniform_joint[i,j] = data[i][2]
        p_MAP_values_Uniform_joint[i,j] = data[i][3]
        q_MAP_values_Uniform_joint[i,j] = data[i][4]
        sigma_MAP_values_Uniform_joint[i,j] = data[i][5]
        tau_MAP_values_Uniform_joint[i,j] = data[i][6]

        ISI_values_Uniform_joint[i,j] = ISI[i]
        computation_time_values_Uniform_joint[i,j] = data[i][8]

        entropy_values_Uniform_joint[i,j] = data[i][9]
    end
end

tau_entropy_Uniform_joint = mean(tau_entropy_values_Uniform_joint,dims=2)
entropy_Uniform_joint = mean(entropy_values_Uniform_joint,dims=2)

N_MAP_Uniform_joint = mean(N_MAP_values_Uniform_joint,dims=2)
p_MAP_Uniform_joint = mean(p_MAP_values_Uniform_joint,dims=2)
q_MAP_Uniform_joint = mean(q_MAP_values_Uniform_joint,dims=2)
sigma_MAP_Uniform_joint = mean(sigma_MAP_values_Uniform_joint,dims=2)
tau_MAP_Uniform_joint = mean(tau_MAP_values_Uniform_joint,dims=2)

# Myopic protocol #########################################################
tau_entropy_values_Myopic = zeros(T,M)
entropy_values_Myopic = zeros(T,M)

N_MAP_values_Myopic = zeros(T,M)
p_MAP_values_Myopic = zeros(T,M)
q_MAP_values_Myopic = zeros(T,M)
sigma_MAP_values_Myopic = zeros(T,M)
tau_MAP_values_Myopic = zeros(T,M)

computation_time_values_Myopic = zeros(T,M)
ISI_values_Myopic = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\Myopic_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][7])
    for i in 1:T
        tau_entropy_values_Myopic[i,j] = data[i][1]

        N_MAP_values_Myopic[i,j] = data[i][2]
        p_MAP_values_Myopic[i,j] = data[i][3]
        q_MAP_values_Myopic[i,j] = data[i][4]
        sigma_MAP_values_Myopic[i,j] = data[i][5]
        tau_MAP_values_Myopic[i,j] = data[i][6]

        ISI_values_Myopic[i,j] = ISI[i]
        computation_time_values_Myopic[i,j] = data[i][8]

        entropy_values_Myopic[i,j] = data[i][9]
    end
end

tau_entropy_Myopic = mean(tau_entropy_values_Myopic,dims=2)
entropy_Myopic = mean(entropy_values_Myopic,dims=2)

N_MAP_Myopic = mean(N_MAP_values_Myopic,dims=2)
p_MAP_Myopic = mean(p_MAP_values_Myopic,dims=2)
q_MAP_Myopic = mean(q_MAP_values_Myopic,dims=2)
sigma_MAP_Myopic = mean(sigma_MAP_values_Myopic,dims=2)
tau_MAP_Myopic = mean(tau_MAP_values_Myopic,dims=2)

###############################################################################
# @save "data.jld"
# @load "data.jld"
N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

###############################################################################

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

savefig("joint_entropy_sim1.png")

################################################################################

RMSE_constant_joint = zeros(T)
RMSE_exponential_joint = zeros(T)
RMSE_Uniform_joint = zeros(T)
RMSE_Myopic = zeros(T)

for i in 1:T

    dist = (1/5)*(
    ((N_MAP_values_constant_joint[i,:] .- N) ./ (maximum(N_MAP_values_constant_joint[i,:])-minimum(N_MAP_values_constant_joint[i,:]))).^2 +
    ((p_MAP_values_constant_joint[i,:] .- p) ./ (maximum(p_MAP_values_constant_joint[i,:])-minimum(p_MAP_values_constant_joint[i,:]))).^2 +
    ((q_MAP_values_constant_joint[i,:] .- q) ./ (maximum(q_MAP_values_constant_joint[i,:])-minimum(q_MAP_values_constant_joint[i,:]))).^2 +
    ((sigma_MAP_values_constant_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_constant_joint[i,:])-minimum(sigma_MAP_values_constant_joint[i,:]))).^2 +
    ((tau_MAP_values_constant_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_constant_joint[i,:])-minimum(tau_MAP_values_constant_joint[i,:]))).^2
    )
    RMSE_constant_joint[i] = sqrt(sum(dist)/M)


    dist = (1/5)*(
    ((N_MAP_values_exponential_joint[i,:] .- N) ./ (maximum(N_MAP_values_exponential_joint[i,:])-minimum(N_MAP_values_exponential_joint[i,:]))).^2 +
    ((p_MAP_values_exponential_joint[i,:] .- p) ./ (maximum(p_MAP_values_exponential_joint[i,:])-minimum(p_MAP_values_exponential_joint[i,:]))).^2 +
    ((q_MAP_values_exponential_joint[i,:] .- q) ./ (maximum(q_MAP_values_exponential_joint[i,:])-minimum(q_MAP_values_exponential_joint[i,:]))).^2 +
    ((sigma_MAP_values_exponential_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_exponential_joint[i,:])-minimum(sigma_MAP_values_exponential_joint[i,:]))).^2 +
    ((tau_MAP_values_exponential_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_exponential_joint[i,:])-minimum(tau_MAP_values_exponential_joint[i,:]))).^2
    )
    RMSE_exponential_joint[i] = sqrt(sum(dist)/M)


    dist = (1/5)*(
    ((N_MAP_values_Uniform_joint[i,:] .- N) ./ (maximum(N_MAP_values_Uniform_joint[i,:])-minimum(N_MAP_values_Uniform_joint[i,:]))).^2 +
    ((p_MAP_values_Uniform_joint[i,:] .- p) ./ (maximum(p_MAP_values_Uniform_joint[i,:])-minimum(p_MAP_values_Uniform_joint[i,:]))).^2 +
    ((q_MAP_values_Uniform_joint[i,:] .- q) ./ (maximum(q_MAP_values_Uniform_joint[i,:])-minimum(q_MAP_values_Uniform_joint[i,:]))).^2 +
    ((sigma_MAP_values_Uniform_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Uniform_joint[i,:])-minimum(sigma_MAP_values_Uniform_joint[i,:]))).^2 +
    ((tau_MAP_values_Uniform_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_Uniform_joint[i,:])-minimum(tau_MAP_values_Uniform_joint[i,:]))).^2
    )
    RMSE_Uniform_joint[i] = sqrt(sum(dist)/M)

    dist = (1/5)*(
    ((N_MAP_values_Myopic[i,:] .- N) ./ (maximum(N_MAP_values_Myopic[i,:])-minimum(N_MAP_values_Myopic[i,:]))).^2 +
    ((p_MAP_values_Myopic[i,:] .- p) ./ (maximum(p_MAP_values_Myopic[i,:])-minimum(p_MAP_values_Myopic[i,:]))).^2 +
    ((q_MAP_values_Myopic[i,:] .- q) ./ (maximum(q_MAP_values_Myopic[i,:])-minimum(q_MAP_values_Myopic[i,:]))).^2 +
    ((sigma_MAP_values_Myopic[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Myopic[i,:])-minimum(sigma_MAP_values_Myopic[i,:]))).^2 +
    ((tau_MAP_values_Myopic[i,:] .- tau) ./ (maximum(tau_MAP_values_Myopic[i,:])-minimum(tau_MAP_values_Myopic[i,:]))).^2
    )
    RMSE_Myopic[i] = sqrt(sum(dist)/M)

end

p1 = plot(RMSE_constant_joint,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300)
plot!(RMSE_exponential_joint,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Uniform_joint,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)

savefig("RMSE_joint_sim1.png")

################################################################################
delta_time_myopic = []

count = 0
for j in 1:100
    for i in 3:200
        append!(delta_time_myopic,ISI_values_Myopic[i,j]-computation_time_values_Myopic[i,j])
    end
end

count/length(delta_time_myopic)

# l = @layout [a ; b]
p1 = histogram(delta_time_myopic,label=false,xlabel="ISI - computation time [s]",
ylabel="# of stimulations")
vline!([0],label=nothing,linewidth = 2,linestyle=:dot,color="red")
p=plot(p1, size = (400, 170),dpi=300)
savefig("computation_time_joint_sim1.png")

################################################################################

