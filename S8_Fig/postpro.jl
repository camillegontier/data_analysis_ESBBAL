using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 200 #Number of data points
M = 400 #Number of simulations

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
        computation_time_values_Myopic[i,j] = data[i][7].time

        entropy_values_Myopic[i,j] = data[i][8]
    end
end

entropy_Myopic = mean(entropy_values_Myopic,dims=2)

N_MAP_Myopic = mean(N_MAP_values_Myopic,dims=2)
p_MAP_Myopic = mean(p_MAP_values_Myopic,dims=2)
q_MAP_Myopic = mean(q_MAP_values_Myopic,dims=2)
sigma_MAP_Myopic = mean(sigma_MAP_values_Myopic,dims=2)
tau_MAP_Myopic = mean(tau_MAP_values_Myopic,dims=2)

# Myopic strat #########################################################
entropy_values_Myopic_strat = zeros(T,M)

N_MAP_values_Myopic_strat = zeros(T,M)
p_MAP_values_Myopic_strat = zeros(T,M)
q_MAP_values_Myopic_strat = zeros(T,M)
sigma_MAP_values_Myopic_strat = zeros(T,M)
tau_MAP_values_Myopic_strat = zeros(T,M)

computation_time_values_Myopic_strat = zeros(T,M)
ISI_values_Myopic_strat = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\Myopic_strat_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Myopic_strat[i,j] = data[i][1]
        p_MAP_values_Myopic_strat[i,j] = data[i][2]
        q_MAP_values_Myopic_strat[i,j] = data[i][3]
        sigma_MAP_values_Myopic_strat[i,j] = data[i][4]
        tau_MAP_values_Myopic_strat[i,j] = data[i][5]

        ISI_values_Myopic_strat[i,j] = ISI[i]
        computation_time_values_Myopic_strat[i,j] = data[i][7].time

        entropy_values_Myopic_strat[i,j] = data[i][8]
    end
end

entropy_Myopic_strat = mean(entropy_values_Myopic_strat,dims=2)

N_MAP_Myopic_strat = mean(N_MAP_values_Myopic_strat,dims=2)
p_MAP_Myopic_strat = mean(p_MAP_values_Myopic_strat,dims=2)
q_MAP_Myopic_strat = mean(q_MAP_values_Myopic_strat,dims=2)
sigma_MAP_Myopic_strat = mean(sigma_MAP_values_Myopic_strat,dims=2)
tau_MAP_Myopic_strat = mean(tau_MAP_values_Myopic_strat,dims=2)

# Constant protocol #########################################################

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
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_constant_joint[i,j] = data[i][1]
        p_MAP_values_constant_joint[i,j] = data[i][2]
        q_MAP_values_constant_joint[i,j] = data[i][3]
        sigma_MAP_values_constant_joint[i,j] = data[i][4]
        tau_MAP_values_constant_joint[i,j] = data[i][5]

        ISI_values_constant_joint[i,j] = ISI[i]
        computation_time_values_constant_joint[i,j] = data[i][7].time

        entropy_values_constant_joint[i,j] = data[i][8]
    end
end

entropy_constant_joint = mean(entropy_values_constant_joint,dims=2)

N_MAP_constant_joint = mean(N_MAP_values_constant_joint,dims=2)
p_MAP_constant_joint = mean(p_MAP_values_constant_joint,dims=2)
q_MAP_constant_joint = mean(q_MAP_values_constant_joint,dims=2)
sigma_MAP_constant_joint = mean(sigma_MAP_values_constant_joint,dims=2)
tau_MAP_constant_joint = mean(tau_MAP_values_constant_joint,dims=2)

# Exponential protocol #########################################################

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
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_exponential_joint[i,j] = data[i][1]
        p_MAP_values_exponential_joint[i,j] = data[i][2]
        q_MAP_values_exponential_joint[i,j] = data[i][3]
        sigma_MAP_values_exponential_joint[i,j] = data[i][4]
        tau_MAP_values_exponential_joint[i,j] = data[i][5]

        ISI_values_exponential_joint[i,j] = ISI[i]
        computation_time_values_exponential_joint[i,j] = data[i][7].time

        entropy_values_exponential_joint[i,j] = data[i][8]
    end
end

entropy_exponential_joint = mean(entropy_values_exponential_joint,dims=2)

N_MAP_exponential_joint = mean(N_MAP_values_exponential_joint,dims=2)
p_MAP_exponential_joint = mean(p_MAP_values_exponential_joint,dims=2)
q_MAP_exponential_joint = mean(q_MAP_values_exponential_joint,dims=2)
sigma_MAP_exponential_joint = mean(sigma_MAP_values_exponential_joint,dims=2)
tau_MAP_exponential_joint = mean(tau_MAP_values_exponential_joint,dims=2)

# Uniform protocol #########################################################
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
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Uniform_joint[i,j] = data[i][1]
        p_MAP_values_Uniform_joint[i,j] = data[i][2]
        q_MAP_values_Uniform_joint[i,j] = data[i][3]
        sigma_MAP_values_Uniform_joint[i,j] = data[i][4]
        tau_MAP_values_Uniform_joint[i,j] = data[i][5]

        ISI_values_Uniform_joint[i,j] = ISI[i]
        computation_time_values_Uniform_joint[i,j] = data[i][7].time

        entropy_values_Uniform_joint[i,j] = data[i][8]
    end
end

entropy_Uniform_joint = mean(entropy_values_Uniform_joint,dims=2)

N_MAP_Uniform_joint = mean(N_MAP_values_Uniform_joint,dims=2)
p_MAP_Uniform_joint = mean(p_MAP_values_Uniform_joint,dims=2)
q_MAP_Uniform_joint = mean(q_MAP_values_Uniform_joint,dims=2)
sigma_MAP_Uniform_joint = mean(sigma_MAP_values_Uniform_joint,dims=2)
tau_MAP_Uniform_joint = mean(tau_MAP_values_Uniform_joint,dims=2)

# Constant strat protocol #########################################################

entropy_values_constant_strat = zeros(T,M)

N_MAP_values_constant_strat = zeros(T,M)
p_MAP_values_constant_strat = zeros(T,M)
q_MAP_values_constant_strat = zeros(T,M)
sigma_MAP_values_constant_strat = zeros(T,M)
tau_MAP_values_constant_strat = zeros(T,M)

computation_time_values_constant_strat = zeros(T,M)
ISI_values_constant_strat = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\constant_strat_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_constant_strat[i,j] = data[i][1]
        p_MAP_values_constant_strat[i,j] = data[i][2]
        q_MAP_values_constant_strat[i,j] = data[i][3]
        sigma_MAP_values_constant_strat[i,j] = data[i][4]
        tau_MAP_values_constant_strat[i,j] = data[i][5]

        ISI_values_constant_strat[i,j] = ISI[i]
        computation_time_values_constant_strat[i,j] = data[i][7].time

        entropy_values_constant_strat[i,j] = data[i][8]
    end
end

entropy_constant_strat = mean(entropy_values_constant_strat,dims=2)

N_MAP_constant_strat = mean(N_MAP_values_constant_strat,dims=2)
p_MAP_constant_strat = mean(p_MAP_values_constant_strat,dims=2)
q_MAP_constant_strat = mean(q_MAP_values_constant_strat,dims=2)
sigma_MAP_constant_strat = mean(sigma_MAP_values_constant_strat,dims=2)
tau_MAP_constant_strat = mean(tau_MAP_values_constant_strat,dims=2)

# Exponential strat protocol #########################################################

entropy_values_exponential_strat = zeros(T,M)

N_MAP_values_exponential_strat = zeros(T,M)
p_MAP_values_exponential_strat = zeros(T,M)
q_MAP_values_exponential_strat = zeros(T,M)
sigma_MAP_values_exponential_strat = zeros(T,M)
tau_MAP_values_exponential_strat = zeros(T,M)

computation_time_values_exponential_strat = zeros(T,M)
ISI_values_exponential_strat = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\exponential_strat_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_exponential_strat[i,j] = data[i][1]
        p_MAP_values_exponential_strat[i,j] = data[i][2]
        q_MAP_values_exponential_strat[i,j] = data[i][3]
        sigma_MAP_values_exponential_strat[i,j] = data[i][4]
        tau_MAP_values_exponential_strat[i,j] = data[i][5]

        ISI_values_exponential_strat[i,j] = ISI[i]
        computation_time_values_exponential_strat[i,j] = data[i][7].time

        entropy_values_exponential_strat[i,j] = data[i][8]
    end
end

entropy_exponential_strat = mean(entropy_values_exponential_strat,dims=2)

N_MAP_exponential_strat = mean(N_MAP_values_exponential_strat,dims=2)
p_MAP_exponential_strat = mean(p_MAP_values_exponential_strat,dims=2)
q_MAP_exponential_strat = mean(q_MAP_values_exponential_strat,dims=2)
sigma_MAP_exponential_strat = mean(sigma_MAP_values_exponential_strat,dims=2)
tau_MAP_exponential_strat = mean(tau_MAP_values_exponential_strat,dims=2)

# Uniform strat protocol #########################################################
entropy_values_Uniform_strat = zeros(T,M)

N_MAP_values_Uniform_strat = zeros(T,M)
p_MAP_values_Uniform_strat = zeros(T,M)
q_MAP_values_Uniform_strat = zeros(T,M)
sigma_MAP_values_Uniform_strat = zeros(T,M)
tau_MAP_values_Uniform_strat = zeros(T,M)

computation_time_values_Uniform_strat = zeros(T,M)
ISI_values_Uniform_strat = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res\\Uniform_strat_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Uniform_strat[i,j] = data[i][1]
        p_MAP_values_Uniform_strat[i,j] = data[i][2]
        q_MAP_values_Uniform_strat[i,j] = data[i][3]
        sigma_MAP_values_Uniform_strat[i,j] = data[i][4]
        tau_MAP_values_Uniform_strat[i,j] = data[i][5]

        ISI_values_Uniform_strat[i,j] = ISI[i]
        computation_time_values_Uniform_strat[i,j] = data[i][7].time

        entropy_values_Uniform_strat[i,j] = data[i][8]
    end
end

entropy_Uniform_strat = mean(entropy_values_Uniform_strat,dims=2)

N_MAP_Uniform_strat = mean(N_MAP_values_Uniform_strat,dims=2)
p_MAP_Uniform_strat = mean(p_MAP_values_Uniform_strat,dims=2)
q_MAP_Uniform_strat = mean(q_MAP_values_Uniform_strat,dims=2)
sigma_MAP_Uniform_strat = mean(sigma_MAP_values_Uniform_strat,dims=2)
tau_MAP_Uniform_strat = mean(tau_MAP_values_Uniform_strat,dims=2)


###############################################################################

N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

@save "data.jld"
@load "data.jld"

###############################################################################

p6 = plot(entropy_constant_joint[1:T],label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300,legend=:topright,
ribbon=std(entropy_values_constant_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Uniform_joint[1:T],label="Best Uniform",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_Uniform_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_exponential_joint[1:T],label="Best Exponential",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_exponential_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic[1:T],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,
color="black",line = :dash,
ribbon=std(entropy_values_Myopic[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

#savefig("exact_entropy.png")

p6 = plot(entropy_constant_strat[1:T],label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300,legend=:topright,
ribbon=std(entropy_values_constant_strat[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Uniform_strat[1:T],label="Best Uniform",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_Uniform_strat[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_exponential_strat[1:T],label="Best Exponential",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_exponential_strat[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic_strat[1:T],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,
color="black",line = :dash,
ribbon=std(entropy_values_Myopic_strat[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

##########
p7 = plot(entropy_Myopic[1:T],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,
color="black",line = :dash,
ribbon=std(entropy_values_Myopic[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic_strat[1:T],label="strat",linewidth = 2, thickness_scaling = 1,
color="grey",line = :dash,
ribbon=std(entropy_values_Myopic_strat[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

##########
p1 = plot(entropy_constant_joint,label="Constant (Multinomial)",linewidth = 2, thickness_scaling = 1,ylabel="Posterior entropy [bit]",
size = (400, 170), dpi=300,
color=palette(:default)[1])
plot!(entropy_constant_strat,label="Constant (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[1],alpha = 0.5,legend=:topright)

p2 = plot(entropy_Uniform_joint,label="Uniform (Multinomial)",linewidth = 2, thickness_scaling = 1,
size = (400, 170), dpi=300,
color=palette(:default)[2])
plot!(entropy_Uniform_strat,label="Uniform (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[2],alpha = 0.5,legend=:topright)

p3 = plot(entropy_exponential_joint,label="Exponential (Multinomial)",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior entropy [bit]",
size = (400, 170), dpi=300,
color=palette(:default)[3])
plot!(entropy_exponential_strat,label="Exponential (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[3],alpha = 0.5)

p4 = plot(entropy_Myopic,label="ESB-BAL (Multinomial)",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",
size = (400, 170), dpi=300,
color="black",line = :dash)
plot!(entropy_Myopic_strat,label="ESB-BAL (Stratified)",linewidth = 2, thickness_scaling = 1,
color="black",alpha = 0.5,line = :dash)

display(plot(p1,p2, p3, p4, layout = 4,dpi=300,size=(600,450)))
savefig("entropy_stratified.png")


#################################################################################


RMSE_constant_joint = zeros(T)
RMSE_exponential_joint = zeros(T)
RMSE_Uniform_joint = zeros(T)
RMSE_Myopic = zeros(T)

RMSE_constant_strat = zeros(T)
RMSE_exponential_strat = zeros(T)
RMSE_Uniform_strat = zeros(T)
RMSE_Myopic_strat = zeros(T)

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_constant_joint[i,:] .- N) ./ (maximum(N_MAP_values_constant_joint[i,:])-minimum(N_MAP_values_constant_joint[i,:]))).^2 +
    ((p_MAP_values_constant_joint[i,:] .- p) ./ (maximum(p_MAP_values_constant_joint[i,:])-minimum(p_MAP_values_constant_joint[i,:]))).^2 +
    ((q_MAP_values_constant_joint[i,:] .- q) ./ (maximum(q_MAP_values_constant_joint[i,:])-minimum(q_MAP_values_constant_joint[i,:]))).^2 +
    ((sigma_MAP_values_constant_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_constant_joint[i,:])-minimum(sigma_MAP_values_constant_joint[i,:]))).^2 +
    ((tau_MAP_values_constant_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_constant_joint[i,:])-minimum(tau_MAP_values_constant_joint[i,:]))).^2
    )
    RMSE_constant_joint[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_exponential_joint[i,:] .- N) ./ (maximum(N_MAP_values_exponential_joint[i,:])-minimum(N_MAP_values_exponential_joint[i,:]))).^2 +
    ((p_MAP_values_exponential_joint[i,:] .- p) ./ (maximum(p_MAP_values_exponential_joint[i,:])-minimum(p_MAP_values_exponential_joint[i,:]))).^2 +
    ((q_MAP_values_exponential_joint[i,:] .- q) ./ (maximum(q_MAP_values_exponential_joint[i,:])-minimum(q_MAP_values_exponential_joint[i,:]))).^2 +
    ((sigma_MAP_values_exponential_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_exponential_joint[i,:])-minimum(sigma_MAP_values_exponential_joint[i,:]))).^2 +
    ((tau_MAP_values_exponential_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_exponential_joint[i,:])-minimum(tau_MAP_values_exponential_joint[i,:]))).^2
    )
    RMSE_exponential_joint[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Uniform_joint[i,:] .- N) ./ (maximum(N_MAP_values_Uniform_joint[i,:])-minimum(N_MAP_values_Uniform_joint[i,:]))).^2 +
    ((p_MAP_values_Uniform_joint[i,:] .- p) ./ (maximum(p_MAP_values_Uniform_joint[i,:])-minimum(p_MAP_values_Uniform_joint[i,:]))).^2 +
    ((q_MAP_values_Uniform_joint[i,:] .- q) ./ (maximum(q_MAP_values_Uniform_joint[i,:])-minimum(q_MAP_values_Uniform_joint[i,:]))).^2 +
    ((sigma_MAP_values_Uniform_joint[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Uniform_joint[i,:])-minimum(sigma_MAP_values_Uniform_joint[i,:]))).^2 +
    ((tau_MAP_values_Uniform_joint[i,:] .- tau) ./ (maximum(tau_MAP_values_Uniform_joint[i,:])-minimum(tau_MAP_values_Uniform_joint[i,:]))).^2
    )
    RMSE_Uniform_joint[i] = sqrt(sum(dist)/M)
end
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
    ((N_MAP_values_constant_strat[i,:] .- N) ./ (maximum(N_MAP_values_constant_strat[i,:])-minimum(N_MAP_values_constant_strat[i,:]))).^2 +
    ((p_MAP_values_constant_strat[i,:] .- p) ./ (maximum(p_MAP_values_constant_strat[i,:])-minimum(p_MAP_values_constant_strat[i,:]))).^2 +
    ((q_MAP_values_constant_strat[i,:] .- q) ./ (maximum(q_MAP_values_constant_strat[i,:])-minimum(q_MAP_values_constant_strat[i,:]))).^2 +
    ((sigma_MAP_values_constant_strat[i,:] .- sigma) ./ (maximum(sigma_MAP_values_constant_strat[i,:])-minimum(sigma_MAP_values_constant_strat[i,:]))).^2 +
    ((tau_MAP_values_constant_strat[i,:] .- tau) ./ (maximum(tau_MAP_values_constant_strat[i,:])-minimum(tau_MAP_values_constant_strat[i,:]))).^2
    )
    RMSE_constant_strat[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_exponential_strat[i,:] .- N) ./ (maximum(N_MAP_values_exponential_strat[i,:])-minimum(N_MAP_values_exponential_strat[i,:]))).^2 +
    ((p_MAP_values_exponential_strat[i,:] .- p) ./ (maximum(p_MAP_values_exponential_strat[i,:])-minimum(p_MAP_values_exponential_strat[i,:]))).^2 +
    ((q_MAP_values_exponential_strat[i,:] .- q) ./ (maximum(q_MAP_values_exponential_strat[i,:])-minimum(q_MAP_values_exponential_strat[i,:]))).^2 +
    ((sigma_MAP_values_exponential_strat[i,:] .- sigma) ./ (maximum(sigma_MAP_values_exponential_strat[i,:])-minimum(sigma_MAP_values_exponential_strat[i,:]))).^2 +
    ((tau_MAP_values_exponential_strat[i,:] .- tau) ./ (maximum(tau_MAP_values_exponential_strat[i,:])-minimum(tau_MAP_values_exponential_strat[i,:]))).^2
    )
    RMSE_exponential_strat[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Uniform_strat[i,:] .- N) ./ (maximum(N_MAP_values_Uniform_strat[i,:])-minimum(N_MAP_values_Uniform_strat[i,:]))).^2 +
    ((p_MAP_values_Uniform_strat[i,:] .- p) ./ (maximum(p_MAP_values_Uniform_strat[i,:])-minimum(p_MAP_values_Uniform_strat[i,:]))).^2 +
    ((q_MAP_values_Uniform_strat[i,:] .- q) ./ (maximum(q_MAP_values_Uniform_strat[i,:])-minimum(q_MAP_values_Uniform_strat[i,:]))).^2 +
    ((sigma_MAP_values_Uniform_strat[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Uniform_strat[i,:])-minimum(sigma_MAP_values_Uniform_strat[i,:]))).^2 +
    ((tau_MAP_values_Uniform_strat[i,:] .- tau) ./ (maximum(tau_MAP_values_Uniform_strat[i,:])-minimum(tau_MAP_values_Uniform_strat[i,:]))).^2
    )
    RMSE_Uniform_strat[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Myopic_strat[i,:] .- N) ./ (maximum(N_MAP_values_Myopic_strat[i,:])-minimum(N_MAP_values_Myopic_strat[i,:]))).^2 +
    ((p_MAP_values_Myopic_strat[i,:] .- p) ./ (maximum(p_MAP_values_Myopic_strat[i,:])-minimum(p_MAP_values_Myopic_strat[i,:]))).^2 +
    ((q_MAP_values_Myopic_strat[i,:] .- q) ./ (maximum(q_MAP_values_Myopic_strat[i,:])-minimum(q_MAP_values_Myopic_strat[i,:]))).^2 +
    ((sigma_MAP_values_Myopic_strat[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Myopic_strat[i,:])-minimum(sigma_MAP_values_Myopic_strat[i,:]))).^2 +
    ((tau_MAP_values_Myopic_strat[i,:] .- tau) ./ (maximum(tau_MAP_values_Myopic_strat[i,:])-minimum(tau_MAP_values_Myopic_strat[i,:]))).^2
    )
    RMSE_Myopic_strat[i] = sqrt(sum(dist)/M)
end


p1 = plot(RMSE_constant_joint,label="Constant (Multinomial)",linewidth = 2, thickness_scaling = 1,ylabel="RMSE",
size = (400, 170), dpi=300,
color=palette(:default)[1])
plot!(RMSE_constant_strat,label="Constant (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[1],alpha = 0.5)

p2 = plot(RMSE_Uniform_joint,label="Uniform (Multinomial)",linewidth = 2, thickness_scaling = 1,
size = (400, 170), dpi=300,
color=palette(:default)[2])
plot!(RMSE_Uniform_strat,label="Uniform (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[2],alpha = 0.5)

p3 = plot(RMSE_exponential_joint,label="Exponential (Multinomial)",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300,
color=palette(:default)[3])
plot!(RMSE_exponential_strat,label="Exponential (Stratified)",linewidth = 2, thickness_scaling = 1,
color=palette(:default)[3],alpha = 0.5)

p4 = plot(RMSE_Myopic,label="ESB-BAL (Multinomial)",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",
size = (400, 170), dpi=300,
color="black",line = :dash)
plot!(RMSE_Myopic_strat,label="ESB-BAL (Stratified)",linewidth = 2, thickness_scaling = 1,
color="black",alpha = 0.5,line = :dash)

display(plot(p1,p2, p3, p4, layout = 4,dpi=300,size=(600,450)))
savefig("RMSE_stratified.png")


##########################################################
##########################################################
##########################################################
##########################################################
##########################################################
##########################################################


p1 = plot(RMSE_constant_joint,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300, legend = false)
plot!(RMSE_Uniform_joint,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_exponential_joint,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)
#plot!(RMSE_fixed_optimal,label="Fixed optimal",linewidth = 2, thickness_scaling = 1,color="purple",line = :dash)

p1 = plot(RMSE_constant_strat,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300, legend = false)
plot!(RMSE_Uniform_strat,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_exponential_strat,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic_strat,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)
#plot!(RMSE_fixed_optimal,label="Fixed optimal",linewidth = 2, thickness_scaling = 1,color="purple",line = :dash)


#################################################################################


RMSE_constant_joint = zeros(T)
RMSE_exponential_joint = zeros(T)
RMSE_Uniform_joint = zeros(T)
RMSE_Myopic = zeros(T)

RMSE_constant_strat = zeros(T)
RMSE_exponential_strat = zeros(T)
RMSE_Uniform_strat = zeros(T)
RMSE_Myopic_strat = zeros(T)

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_constant_joint[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_constant_joint[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_constant_joint[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_constant_joint[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_constant_joint[i,:] .- tau) ./ tau).^2
    )
    RMSE_constant_joint[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_exponential_joint[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_exponential_joint[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_exponential_joint[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_exponential_joint[i,:] .- sigma) ./sigma).^2 +
    ((tau_MAP_values_exponential_joint[i,:] .- tau) ./ tau).^2
    )
    RMSE_exponential_joint[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Uniform_joint[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_Uniform_joint[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_Uniform_joint[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_Uniform_joint[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_Uniform_joint[i,:] .- tau) ./ tau).^2
    )
    RMSE_Uniform_joint[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Myopic[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_Myopic[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_Myopic[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_Myopic[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_Myopic[i,:] .- tau) ./ tau).^2
    )
    RMSE_Myopic[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_constant_strat[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_constant_strat[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_constant_strat[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_constant_strat[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_constant_strat[i,:] .- tau) ./ tau).^2
    )
    RMSE_constant_strat[i] = sqrt(sum(dist)/M)
end

for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_exponential_strat[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_exponential_strat[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_exponential_strat[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_exponential_strat[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_exponential_strat[i,:] .- tau) ./ tau).^2
    )
    RMSE_exponential_strat[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Uniform_strat[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_Uniform_strat[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_Uniform_strat[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_Uniform_strat[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_Uniform_strat[i,:] .- tau) ./ tau).^2
    )
    RMSE_Uniform_strat[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_Myopic_strat[i,:] .- N) ./ N).^2 +
    ((p_MAP_values_Myopic_strat[i,:] .- p) ./ p).^2 +
    ((q_MAP_values_Myopic_strat[i,:] .- q) ./ q).^2 +
    ((sigma_MAP_values_Myopic_strat[i,:] .- sigma) ./ sigma).^2 +
    ((tau_MAP_values_Myopic_strat[i,:] .- tau) ./ tau).^2
    )
    RMSE_Myopic_strat[i] = sqrt(sum(dist)/M)
end



p1 = plot(RMSE_constant_joint,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300, legend = false)
plot!(RMSE_Uniform_joint,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_exponential_joint,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)
#plot!(RMSE_fixed_optimal,label="Fixed optimal",linewidth = 2, thickness_scaling = 1,color="purple",line = :dash)

p1 = plot(RMSE_constant_strat,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300, legend = false)
plot!(RMSE_Uniform_strat,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_exponential_strat,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic_strat,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)
#plot!(RMSE_fixed_optimal,label="Fixed optimal",linewidth = 2, thickness_scaling = 1,color="purple",line = :dash)
