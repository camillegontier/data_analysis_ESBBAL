using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 200 #Number of data points
M = 400 #Number of simulations

# short protocol #########################################################
entropy_values_short = zeros(T,M)

N_MAP_values_short = zeros(T,M)
p_MAP_values_short = zeros(T,M)
q_MAP_values_short = zeros(T,M)
sigma_MAP_values_short = zeros(T,M)
tau_MAP_values_short = zeros(T,M)

computation_time_values_short = zeros(T,M)
ISI_values_short = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\deterministic_short_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_short[i,j] = data[i][1]
        p_MAP_values_short[i,j] = data[i][2]
        q_MAP_values_short[i,j] = data[i][3]
        sigma_MAP_values_short[i,j] = data[i][4]
        tau_MAP_values_short[i,j] = data[i][5]

        ISI_values_short[i,j] = ISI[i]
        computation_time_values_short[i,j] = data[i][7]

        entropy_values_short[i,j] = data[i][8]
    end
end

entropy_short = mean(entropy_values_short,dims=2)

N_MAP_short = mean(N_MAP_values_short,dims=2)
p_MAP_short = mean(p_MAP_values_short,dims=2)
q_MAP_short = mean(q_MAP_values_short,dims=2)
sigma_MAP_short = mean(sigma_MAP_values_short,dims=2)
tau_MAP_short = mean(tau_MAP_values_short,dims=2)


# long protocol #########################################################
entropy_values_long = zeros(T,M)

N_MAP_values_long = zeros(T,M)
p_MAP_values_long = zeros(T,M)
q_MAP_values_long = zeros(T,M)
sigma_MAP_values_long = zeros(T,M)
tau_MAP_values_long = zeros(T,M)

computation_time_values_long = zeros(T,M)
ISI_values_long = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\deterministic_long_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_long[i,j] = data[i][1]
        p_MAP_values_long[i,j] = data[i][2]
        q_MAP_values_long[i,j] = data[i][3]
        sigma_MAP_values_long[i,j] = data[i][4]
        tau_MAP_values_long[i,j] = data[i][5]

        ISI_values_long[i,j] = ISI[i]
        computation_time_values_long[i,j] = data[i][7]

        entropy_values_long[i,j] = data[i][8]
    end
end

entropy_long = mean(entropy_values_long,dims=2)

N_MAP_long = mean(N_MAP_values_long,dims=2)
p_MAP_long = mean(p_MAP_values_long,dims=2)
q_MAP_long = mean(q_MAP_values_long,dims=2)
sigma_MAP_long = mean(sigma_MAP_values_long,dims=2)
tau_MAP_long = mean(tau_MAP_values_long,dims=2)

# Myopic protocol #########################################################
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
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Myopic[i,j] = data[i][1]
        p_MAP_values_Myopic[i,j] = data[i][2]
        q_MAP_values_Myopic[i,j] = data[i][3]
        sigma_MAP_values_Myopic[i,j] = data[i][4]
        tau_MAP_values_Myopic[i,j] = data[i][5]

        ISI_values_Myopic[i,j] = ISI[i]
        computation_time_values_Myopic[i,j] = data[i][7]

        entropy_values_Myopic[i,j] = data[i][8]
    end
end

entropy_Myopic = mean(entropy_values_Myopic,dims=2)

N_MAP_Myopic = mean(N_MAP_values_Myopic,dims=2)
p_MAP_Myopic = mean(p_MAP_values_Myopic,dims=2)
q_MAP_Myopic = mean(q_MAP_values_Myopic,dims=2)
sigma_MAP_Myopic = mean(sigma_MAP_values_Myopic,dims=2)
tau_MAP_Myopic = mean(tau_MAP_values_Myopic,dims=2)

# batch protocol #########################################################
entropy_values_batch = zeros(T,M)

N_MAP_values_batch = zeros(T,M)
p_MAP_values_batch = zeros(T,M)
q_MAP_values_batch = zeros(T,M)
sigma_MAP_values_batch = zeros(T,M)
tau_MAP_values_batch = zeros(T,M)

computation_time_values_batch = zeros(T,M)
ISI_values_batch = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\batch1_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_batch[i,j] = data[i][1]
        p_MAP_values_batch[i,j] = data[i][2]
        q_MAP_values_batch[i,j] = data[i][3]
        sigma_MAP_values_batch[i,j] = data[i][4]
        tau_MAP_values_batch[i,j] = data[i][5]

        ISI_values_batch[i,j] = ISI[i]
        computation_time_values_batch[i,j] = data[i][7]

        entropy_values_batch[i,j] = data[i][8]
    end
end

entropy_batch = mean(entropy_values_batch,dims=2)

N_MAP_batch = mean(N_MAP_values_batch,dims=2)
p_MAP_batch = mean(p_MAP_values_batch,dims=2)
q_MAP_batch = mean(q_MAP_values_batch,dims=2)
sigma_MAP_batch = mean(sigma_MAP_values_batch,dims=2)
tau_MAP_batch = mean(tau_MAP_values_batch,dims=2)

###############################################################################

N = 47
p = 0.27
q = 1.0
sigma = 0.5
tau = 0.17

@save "data_batch.jld"
@load "data_batch.jld"

###############################################################################

p6 = plot(entropy_long[1:T],label="Deterministic (long)",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300,legend=:topright,
ribbon=std(entropy_values_long[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_short[1:T],label="Deterministic (short)",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_short[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic[1:T],label="ESB-BAL (myopic)",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_Myopic[1:T,:],dims=2)/sqrt(M),fillalpha=.5,
color="black",line = :dash)

plot!(entropy_batch[1:T],label="ESB-BAL (batch)",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_batch[1:T,:],dims=2)/sqrt(M),fillalpha=.5,
color="black",line = :dot)

savefig("batch_entropy.png")

p6 = plot(entropy_Myopic[1:T],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300,legend=:topright,
ribbon=std(entropy_values_Myopic[1:T,:],dims=2)/sqrt(M),fillalpha=.5,
color="black",line = :dash)

#################################################################################


RMSE_Myopic = zeros(T)
RMSE_batch = zeros(T)
RMSE_short = zeros(T)
RMSE_long = zeros(T)

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Myopic[i,:] .- N) ./ (maximum(N_MAP_values_Myopic[i,:])-minimum(N_MAP_values_Myopic[i,:]))).^2 +
    ((p_MAP_values_Myopic[i,:] .- p) ./ (maximum(p_MAP_values_Myopic[i,:])-minimum(p_MAP_values_Myopic[i,:]))).^2 +
    ((q_MAP_values_Myopic[i,:] .- q) ./ (maximum(q_MAP_values_Myopic[i,:])-minimum(q_MAP_values_Myopic[i,:]))).^2 +
    ((sigma_MAP_values_Myopic[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Myopic[i,:])-minimum(sigma_MAP_values_Myopic[i,:]))).^2 +
    ((tau_MAP_values_Myopic[i,:] .- tau) ./ (maximum(tau_MAP_values_Myopic[i,:])-minimum(tau_MAP_values_Myopic[i,:]))).^2
    )
    RMSE_Myopic[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_short[i,:] .- N) ./ (maximum(N_MAP_values_short[i,:])-minimum(N_MAP_values_short[i,:]))).^2 +
    ((p_MAP_values_short[i,:] .- p) ./ (maximum(p_MAP_values_short[i,:])-minimum(p_MAP_values_short[i,:]))).^2 +
    ((q_MAP_values_short[i,:] .- q) ./ (maximum(q_MAP_values_short[i,:])-minimum(q_MAP_values_short[i,:]))).^2 +
    ((sigma_MAP_values_short[i,:] .- sigma) ./ (maximum(sigma_MAP_values_short[i,:])-minimum(sigma_MAP_values_short[i,:]))).^2 +
    ((tau_MAP_values_short[i,:] .- tau) ./ (maximum(tau_MAP_values_short[i,:])-minimum(tau_MAP_values_short[i,:]))).^2
    )
    RMSE_short[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_long[i,:] .- N) ./ (maximum(N_MAP_values_long[i,:])-minimum(N_MAP_values_long[i,:]))).^2 +
    ((p_MAP_values_long[i,:] .- p) ./ (maximum(p_MAP_values_long[i,:])-minimum(p_MAP_values_long[i,:]))).^2 +
    ((q_MAP_values_long[i,:] .- q) ./ (maximum(q_MAP_values_long[i,:])-minimum(q_MAP_values_long[i,:]))).^2 +
    ((sigma_MAP_values_long[i,:] .- sigma) ./ (maximum(sigma_MAP_values_long[i,:])-minimum(sigma_MAP_values_long[i,:]))).^2 +
    ((tau_MAP_values_long[i,:] .- tau) ./ (maximum(tau_MAP_values_long[i,:])-minimum(tau_MAP_values_long[i,:]))).^2
    )
    RMSE_long[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_batch[i,:] .- N) ./ (maximum(N_MAP_values_batch[i,:])-minimum(N_MAP_values_batch[i,:]))).^2 +
    ((p_MAP_values_batch[i,:] .- p) ./ (maximum(p_MAP_values_batch[i,:])-minimum(p_MAP_values_batch[i,:]))).^2 +
    ((q_MAP_values_batch[i,:] .- q) ./ (maximum(q_MAP_values_batch[i,:])-minimum(q_MAP_values_batch[i,:]))).^2 +
    ((sigma_MAP_values_batch[i,:] .- sigma) ./ (maximum(sigma_MAP_values_batch[i,:])-minimum(sigma_MAP_values_batch[i,:]))).^2 +
    ((tau_MAP_values_batch[i,:] .- tau) ./ (maximum(tau_MAP_values_batch[i,:])-minimum(tau_MAP_values_batch[i,:]))).^2
    )
    RMSE_batch[i] = sqrt(sum(dist)/M)
end



p1 = plot(RMSE_Myopic,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300,color="black",line = :dash)

plot!(RMSE_long,label="Long",linewidth = 2, thickness_scaling = 1,color="grey",line = :dash)
plot!(RMSE_short,label="Short",linewidth = 2, thickness_scaling = 1,color="grey",line = :dot)
plot!(RMSE_batch,label="Batch",linewidth = 2, thickness_scaling = 1,color="grey",line = :dot)

savefig("exact_RMSE.png")
