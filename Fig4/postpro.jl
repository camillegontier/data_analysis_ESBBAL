using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 400 #Number of data points
M = 100 #Number of simulations

# Ground truth parameters
N = 7
p = 0.6
q = 1.0
sigma = 0.2
tau = 0.25

#########################################################
ISI_values_00 = zeros(T,M)
entropy_values_00 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_00_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_00[i,j] = data[i][1]
        ISI_values_00[i,j] = ISI[i]
    end
end
entropy_00 = mean(entropy_values_00,dims=2)

#########################################################
ISI_values_01 = zeros(T,M)
entropy_values_01 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_01_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_01[i,j] = data[i][1]
        ISI_values_01[i,j] = ISI[i]
    end
end
entropy_01 = mean(entropy_values_01,dims=2)

#########################################################
ISI_values_02 = zeros(T,M)
entropy_values_02 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_02_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_02[i,j] = data[i][1]
        ISI_values_02[i,j] = ISI[i]
    end
end
entropy_02 = mean(entropy_values_02,dims=2)

#########################################################
ISI_values_03 = zeros(T,M)
entropy_values_03 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_03_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_03[i,j] = data[i][1]
        ISI_values_03[i,j] = ISI[i]
    end
end
entropy_03 = mean(entropy_values_03,dims=2)
#########################################################
ISI_values_04 = zeros(T,M)
entropy_values_04 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_04_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_04[i,j] = data[i][1]
        ISI_values_04[i,j] = ISI[i]
    end
end
entropy_04 = mean(entropy_values_04,dims=2)

#########################################################
ISI_values_05 = zeros(T,M)
entropy_values_05 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_05_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_05[i,j] = data[i][1]
        ISI_values_05[i,j] = ISI[i]
    end
end
entropy_05 = mean(entropy_values_05,dims=2)
#########################################################
ISI_values_06 = zeros(T,M)
entropy_values_06 = zeros(T,M)
for j in 1:M
    print(j)
    obj = load(string("res3\\main_06_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        entropy_values_06[i,j] = data[i][1]
        ISI_values_06[i,j] = ISI[i]
    end
end
entropy_06 = mean(entropy_values_06,dims=2)

# #########################################################


@save "data3.jld"
# @load "data3.jld"
################################################################################
default(palette = palette(:viridis,8, rev = false))
p2 = plot(entropy_00[1:400],xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",label=L"\eta = 0",linewidth = 2, thickness_scaling = 1,
 size = (400, 300), dpi=300)
plot!(entropy_01[1:400],label=L"\eta = 0.6",linewidth = 2, thickness_scaling = 1)
plot!(entropy_02[1:400],label=L"\eta = 1.2",linewidth = 2, thickness_scaling = 1)
plot!(entropy_03[1:400],label=L"\eta = 1.8",linewidth = 2, thickness_scaling = 1)
plot!(entropy_04[1:400],label=L"\eta = 2.4",linewidth = 2, thickness_scaling = 1)
plot!(entropy_05[1:400],label=L"\eta = 3.0",linewidth = 2, thickness_scaling = 1)
# plot!(entropy_06[1:400],label=L"\eta = 0.6",linewidth = 2, thickness_scaling = 1)

###############################################################################
TT = 200
delta_cum = zeros(T,M,7)
time_values = zeros(TT,M,7)

for j in 1:M

    delta_cum[:,j,1] = cumsum(ISI_values_00[:,j])
    delta_cum[:,j,2] = cumsum(ISI_values_01[:,j])
    delta_cum[:,j,3] = cumsum(ISI_values_02[:,j])
    delta_cum[:,j,4] = cumsum(ISI_values_03[:,j])
    delta_cum[:,j,5] = cumsum(ISI_values_04[:,j])
    delta_cum[:,j,6] = cumsum(ISI_values_05[:,j])
    delta_cum[:,j,7] = cumsum(ISI_values_06[:,j])
end

times = LinRange(0.05,15,TT)
for i in 1:length(times)
    t = times[i]
    for j in 1:M
        idx = findall(x->x <= t,delta_cum[:,j,1])
        if isempty(idx)
            time_values[i,j,1] = entropy_values_00[1,j]
        else
            time_values[i,j,1] = entropy_values_00[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,2])
        if isempty(idx)
            time_values[i,j,2] = entropy_values_01[1,j]
        else
            time_values[i,j,2] = entropy_values_01[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,3])
        if isempty(idx)
            time_values[i,j,3] = entropy_values_02[1,j]
        else
            time_values[i,j,3] = entropy_values_02[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,4])
        if isempty(idx)
            time_values[i,j,4] = entropy_values_03[1,j]
        else
            time_values[i,j,4] = entropy_values_03[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,5])
        if isempty(idx)
            time_values[i,j,5] = entropy_values_04[1,j]
        else
            time_values[i,j,5] = entropy_values_04[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,6])
        if isempty(idx)
            time_values[i,j,6] = entropy_values_05[1,j]
        else
            time_values[i,j,6] = entropy_values_05[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,7])
        if isempty(idx)
            time_values[i,j,7] = entropy_values_06[1,j]
        else
            time_values[i,j,7] = entropy_values_06[idx[end],j]
        end
    end
end

time_mean= mean(time_values,dims=2)


############################################################################

entropy_decrease = zeros(M,7)
final_time = zeros(M,7)
info_rate = zeros(M,7)
t_final = 200

info = zeros(7)
for j in 1:7
    info[j] = (time_mean[1,1,j]-time_mean[end,1,j])/times[end]
end


default(palette = palette(:viridis,8, rev = false))
p2 = plot(times,time_mean[:,1,1],xlabel="Time [s]",ylabel="Posterior Entropy [bit]", legend = false,linewidth = 2, thickness_scaling = 1, size = (400, 300), dpi=300)
plot!(times,time_mean[:,1,2],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,3],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,4],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,5],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,6],linewidth = 2, legend = false)
plot!([0,0.1,0.2,0.3,0.4,0.5,0.6],info, inset = (1, bbox(0.6,0.0,0.4,0.4)), subplot = 2, legend = false,linewidth = 2,xlabel=L"\eta",ylabel="[bit/s]",color="black")




################################################################################



