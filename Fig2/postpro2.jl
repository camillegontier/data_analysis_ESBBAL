using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 200 #Number of data points
M = 100 #Number of simulations

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
    obj = load(string("res3\\constant_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_constant_joint[i,j] = data[i][1]
        p_MAP_values_constant_joint[i,j] = data[i][2]
        q_MAP_values_constant_joint[i,j] = data[i][3]
        sigma_MAP_values_constant_joint[i,j] = data[i][4]
        tau_MAP_values_constant_joint[i,j] = data[i][5]

        ISI_values_constant_joint[i,j] = ISI[i]
        computation_time_values_constant_joint[i,j] = data[i][7]

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
    obj = load(string("res3\\exponential_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_exponential_joint[i,j] = data[i][1]
        p_MAP_values_exponential_joint[i,j] = data[i][2]
        q_MAP_values_exponential_joint[i,j] = data[i][3]
        sigma_MAP_values_exponential_joint[i,j] = data[i][4]
        tau_MAP_values_exponential_joint[i,j] = data[i][5]

        ISI_values_exponential_joint[i,j] = ISI[i]
        computation_time_values_exponential_joint[i,j] = data[i][7]

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
    obj = load(string("res3\\Uniform_joint_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Uniform_joint[i,j] = data[i][1]
        p_MAP_values_Uniform_joint[i,j] = data[i][2]
        q_MAP_values_Uniform_joint[i,j] = data[i][3]
        sigma_MAP_values_Uniform_joint[i,j] = data[i][4]
        tau_MAP_values_Uniform_joint[i,j] = data[i][5]

        ISI_values_Uniform_joint[i,j] = ISI[i]
        computation_time_values_Uniform_joint[i,j] = data[i][7]

        entropy_values_Uniform_joint[i,j] = data[i][8]
    end
end

entropy_Uniform_joint = mean(entropy_values_Uniform_joint,dims=2)

N_MAP_Uniform_joint = mean(N_MAP_values_Uniform_joint,dims=2)
p_MAP_Uniform_joint = mean(p_MAP_values_Uniform_joint,dims=2)
q_MAP_Uniform_joint = mean(q_MAP_values_Uniform_joint,dims=2)
sigma_MAP_Uniform_joint = mean(sigma_MAP_values_Uniform_joint,dims=2)
tau_MAP_Uniform_joint = mean(tau_MAP_values_Uniform_joint,dims=2)

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
    obj = load(string("res3\\Myopic_",j,".jld"))
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

# Myopic exact protocol #########################################################
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
    obj = load(string("res3\\Myopic_exact1_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_Myopic_exact[i,j] = data[i][1]
        p_MAP_values_Myopic_exact[i,j] = data[i][2]
        q_MAP_values_Myopic_exact[i,j] = data[i][3]
        sigma_MAP_values_Myopic_exact[i,j] = data[i][4]
        tau_MAP_values_Myopic_exact[i,j] = data[i][5]

        ISI_values_Myopic_exact[i,j] = ISI[i]
        computation_time_values_Myopic_exact[i,j] = data[i][7]

        entropy_values_Myopic_exact[i,j] = data[i][8]
    end
end

entropy_Myopic_exact = mean(entropy_values_Myopic_exact,dims=2)

N_MAP_Myopic_exact = mean(N_MAP_values_Myopic_exact,dims=2)
p_MAP_Myopic_exact = mean(p_MAP_values_Myopic_exact,dims=2)
q_MAP_Myopic_exact = mean(q_MAP_values_Myopic_exact,dims=2)
sigma_MAP_Myopic_exact = mean(sigma_MAP_values_Myopic_exact,dims=2)
tau_MAP_Myopic_exact = mean(tau_MAP_values_Myopic_exact,dims=2)

# Fixed optimal protocol #########################################################
entropy_values_fixed_optimal = zeros(T,M)

N_MAP_values_fixed_optimal = zeros(T,M)
p_MAP_values_fixed_optimal = zeros(T,M)
q_MAP_values_fixed_optimal = zeros(T,M)
sigma_MAP_values_fixed_optimal = zeros(T,M)
tau_MAP_values_fixed_optimal = zeros(T,M)

computation_time_values_fixed_optimal = zeros(T,M)
ISI_values_fixed_optimal = zeros(T,M)

for j in 1:M
    print(j)
    obj = load(string("res3\\optimal_fixed_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][6])
    for i in 1:T

        N_MAP_values_fixed_optimal[i,j] = data[i][1]
        p_MAP_values_fixed_optimal[i,j] = data[i][2]
        q_MAP_values_fixed_optimal[i,j] = data[i][3]
        sigma_MAP_values_fixed_optimal[i,j] = data[i][4]
        tau_MAP_values_fixed_optimal[i,j] = data[i][5]

        ISI_values_fixed_optimal[i,j] = ISI[i]
        computation_time_values_fixed_optimal[i,j] = data[i][7]

        entropy_values_fixed_optimal[i,j] = data[i][8]
    end
end

entropy_fixed_optimal = mean(entropy_values_fixed_optimal,dims=2)

N_MAP_fixed_optimal = mean(N_MAP_values_fixed_optimal,dims=2)
p_MAP_fixed_optimal = mean(p_MAP_values_fixed_optimal,dims=2)
q_MAP_fixed_optimal = mean(q_MAP_values_fixed_optimal,dims=2)
sigma_MAP_fixed_optimal = mean(sigma_MAP_values_fixed_optimal,dims=2)
tau_MAP_fixed_optimal = mean(tau_MAP_values_fixed_optimal,dims=2)

###############################################################################
N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

@save "data3.jld"
# @load "data.jld"


###############################################################################

p6 = plot(entropy_constant_joint[1:T],label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",
size = (400, 300), dpi=300, title="Joint distribution",
ribbon=std(entropy_values_constant_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Uniform_joint[1:T],label="Best Uniform",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_Uniform_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_exponential_joint[1:T],label="Best Exponential",linewidth = 2, thickness_scaling = 1,
ribbon=std(entropy_values_exponential_joint[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic[1:T],label="ESB-BAL",linewidth = 2, thickness_scaling = 1,
color="black",line = :dash,
ribbon=std(entropy_values_Myopic[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_Myopic_exact[1:T],label="ESB-BAL (\"exact\")",linewidth = 2, thickness_scaling = 1,
color="grey",line = :dash,
ribbon=std(entropy_values_Myopic_exact[1:T,:],dims=2)/sqrt(M),fillalpha=.5)

plot!(entropy_fixed_optimal[1:T],label="Fixed optimal",linewidth = 2, thickness_scaling = 1,
color="grey",line = :dash,
ribbon=std(entropy_values_fixed_optimal[1:T,:],dims=2)/sqrt(M),fillalpha=.5)



#savefig("joint_entropy_sim1.png")

################################################################################


RMSE_constant_joint = zeros(T)
RMSE_exponential_joint = zeros(T)
RMSE_Uniform_joint = zeros(T)
RMSE_Myopic = zeros(T)
RMSE_Myopic_exact = zeros(T)
RMSE_fixed_optimal = zeros(T)

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
    ((N_MAP_values_Myopic_exact[i,:] .- N) ./ (maximum(N_MAP_values_Myopic_exact[i,:])-minimum(N_MAP_values_Myopic_exact[i,:]))).^2 +
    ((p_MAP_values_Myopic_exact[i,:] .- p) ./ (maximum(p_MAP_values_Myopic_exact[i,:])-minimum(p_MAP_values_Myopic_exact[i,:]))).^2 +
    ((q_MAP_values_Myopic_exact[i,:] .- q) ./ (maximum(q_MAP_values_Myopic_exact[i,:])-minimum(q_MAP_values_Myopic_exact[i,:]))).^2 +
    ((sigma_MAP_values_Myopic_exact[i,:] .- sigma) ./ (maximum(sigma_MAP_values_Myopic_exact[i,:])-minimum(sigma_MAP_values_Myopic_exact[i,:]))).^2 +
    ((tau_MAP_values_Myopic_exact[i,:] .- tau) ./ (maximum(tau_MAP_values_Myopic_exact[i,:])-minimum(tau_MAP_values_Myopic_exact[i,:]))).^2
    )
    RMSE_Myopic_exact[i] = sqrt(sum(dist)/M)
end
for i in 1:T
    dist = (1/5)*(
    ((N_MAP_values_fixed_optimal[i,:] .- N) ./ (maximum(N_MAP_values_fixed_optimal[i,:])-minimum(N_MAP_values_fixed_optimal[i,:]))).^2 +
    ((p_MAP_values_fixed_optimal[i,:] .- p) ./ (maximum(p_MAP_values_fixed_optimal[i,:])-minimum(p_MAP_values_fixed_optimal[i,:]))).^2 +
    ((q_MAP_values_fixed_optimal[i,:] .- q) ./ (maximum(q_MAP_values_fixed_optimal[i,:])-minimum(q_MAP_values_fixed_optimal[i,:]))).^2 +
    ((sigma_MAP_values_fixed_optimal[i,:] .- sigma) ./ (maximum(sigma_MAP_values_fixed_optimal[i,:])-minimum(sigma_MAP_values_fixed_optimal[i,:]))).^2 +
    ((tau_MAP_values_fixed_optimal[i,:] .- tau) ./ (maximum(tau_MAP_values_fixed_optimal[i,:])-minimum(tau_MAP_values_fixed_optimal[i,:]))).^2
    )
    RMSE_fixed_optimal[i] = sqrt(sum(dist)/M)
end

p1 = plot(RMSE_constant_joint,label="Best Constant",linewidth = 2, thickness_scaling = 1,xlabel="# of stimulations",ylabel="RMSE",
size = (400, 170), dpi=300, legend = false)
plot!(RMSE_Uniform_joint,label="Best Uniform",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_exponential_joint,label="Best Exponential",linewidth = 2, thickness_scaling = 1)
plot!(RMSE_Myopic,label="ESB-BAL",linewidth = 2, thickness_scaling = 1,color="black",line = :dash)
plot!(RMSE_Myopic_exact,label="ESB-BAL exact",linewidth = 2, thickness_scaling = 1,color="grey",line = :dash)
plot!(RMSE_fixed_optimal,label="Fixed optimal",linewidth = 2, thickness_scaling = 1,color="grey",line = :dash)

#savefig("RMSE_joint_sim1.png")

################################################################################

x = []
y = []

count = 0
for j in 1:100
    for i in 2:200
        append!(x,ISI_values_Myopic[i,j])
    end
end

layout = @layout [a            _
                  b{0.7w,0.7h} c]

default(fillcolor = :blue, markercolor = :white, grid = false, legend = false)
plot(layout = layout, link = :both, size = (400, 200),  right_margin = -35Plots.px,  top_margin = -20Plots.px)
scatter!(x,y, subplot = 2, framestyle = :box,xlabel="ISI [s]",
ylabel="Computation time [s]")
histogram!([x y], subplot = [1 3], orientation = [:v :h], framestyle = :none)

