using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

M = 400
res_file = "results"

# Constant #####################################################################
hyperparameters = LinRange(0.05,1,30)
res_constant = zeros(length(hyperparameters))
res_constant_tau = zeros(length(hyperparameters))

for i in 1:length(hyperparameters)
    print(i)
    res_constant_temp = []
    res_constant_temp_tau = []
    for j in 1:M
        try
            obj = load(string(res_file,"/constant_",i,"_",j,".jld"))
            data = obj["data"]
            append!(res_constant_temp,data[end][2] - data[1][2])
            append!(res_constant_temp_tau,data[end][1] - data[1][1])
        catch e
        end

    end
    res_constant[i] = mean(res_constant_temp)
    res_constant_tau[i] = mean(res_constant_temp_tau)
end
# Exponential #####################################################################
hyperparameters = LinRange(0.05,1,30)
res_exponential = zeros(length(hyperparameters))
res_exponential_tau = zeros(length(hyperparameters))

for i in 1:length(hyperparameters)
    print(i)
    res_exponential_temp = []
    res_exponential_temp_tau = []
    for j in 1:M
        try
            obj = load(string(res_file,"/exponential_",i,"_",j,".jld"))
            data = obj["data"]
            append!(res_exponential_temp,data[end][2] - data[1][2])
            append!(res_exponential_temp_tau,data[end][1] - data[1][1])
        catch e
        end

    end
    res_exponential[i] = mean(res_exponential_temp)
    res_exponential_tau[i] = mean(res_exponential_temp_tau)
end
# Uniform #####################################################################
hyperparameters = LinRange(0.05,2,30)
res_uniform = zeros(length(hyperparameters))
res_uniform_tau = zeros(length(hyperparameters))

for i in 1:length(hyperparameters)
    print(i)
    res_uniform_temp = []
    res_uniform_temp_tau = []
    for j in 1:M
        try
            obj = load(string(res_file,"/Uniform_",i,"_",j,".jld"))
            data = obj["data"]
            append!(res_uniform_temp,data[end][2] - data[1][2])
            append!(res_uniform_temp_tau,data[end][1] - data[1][1])
        catch e
        end

    end
    try
        res_uniform[i] = mean(res_uniform_temp)
        res_uniform_tau[i] = mean(res_uniform_temp_tau)
    catch e
    end
end

#####################################################################
print(res_constant)
print('\n')
print(res_exponential)
print('\n')
print(res_uniform)
print('\n')

print(res_constant_tau)
print('\n')
print(res_exponential_tau)
print('\n')
print(res_uniform_tau)
print('\n')
