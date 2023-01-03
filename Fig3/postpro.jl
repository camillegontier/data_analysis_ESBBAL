using JLD
using Plots
using StatsBase
using Statistics
using LaTeXStrings

T = 200 #Number of data points
M = 200 #Number of simulations

# nu = 0.0 #########################################################

tau_entropy_values_00 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_00_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_00[i,j] = ISI[i]
        tau_entropy_values_00[i,j] = data[i][1]

    end
end

tau_entropy_00 = mean(tau_entropy_values_00,dims=2)

# nu = 0.1 #########################################################

tau_entropy_values_01 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_01_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_01[i,j] = ISI[i]
        tau_entropy_values_01[i,j] = data[i][1]

    end
end

tau_entropy_01 = mean(tau_entropy_values_01,dims=2)

# nu = 0.2 #########################################################

tau_entropy_values_02 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_02_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_02[i,j] = ISI[i]
        tau_entropy_values_02[i,j] = data[i][1]

    end
end

tau_entropy_02 = mean(tau_entropy_values_02,dims=2)

# nu = 0.3 #########################################################

tau_entropy_values_03 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_03_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_03[i,j] = ISI[i]
        tau_entropy_values_03[i,j] = data[i][1]

    end
end

tau_entropy_03 = mean(tau_entropy_values_03,dims=2)

# nu = 0.4 #########################################################

tau_entropy_values_04 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_04_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_04[i,j] = ISI[i]
        tau_entropy_values_04[i,j] = data[i][1]

    end
end

tau_entropy_04 = mean(tau_entropy_values_04,dims=2)

# nu = 0.5 #########################################################

tau_entropy_values_05 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_05_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_05[i,j] = ISI[i]
        tau_entropy_values_05[i,j] = data[i][1]

    end
end

tau_entropy_05 = mean(tau_entropy_values_05,dims=2)

# nu = 0.6 #########################################################

tau_entropy_values_06 = zeros(T,M)

for j in 1:M
    print(j)

    obj = load(string("res\\main_06_",j,".jld"))
    data = obj["data"]
    ISI = diff(data[1][2])
    for i in 1:T
        ISI_values_06[i,j] = ISI[i]
        tau_entropy_values_06[i,j] = data[i][1]

    end
end

tau_entropy_06 = mean(tau_entropy_values_06,dims=2)

@save "data.jld"
# @load "data.jld"

################################################################################
default(palette = palette(:viridis,9, rev = false))
p2 = plot(tau_entropy_00[1:200],xlabel="# of stimulations",ylabel="Posterior Entropy [bit]",label=L"\eta = 0",linewidth = 2, thickness_scaling = 1,
 size = (400, 300), dpi=300)
plot!(tau_entropy_01[1:200],label=L"\eta = 0.1",linewidth = 2, thickness_scaling = 1)
plot!(tau_entropy_02[1:200],label=L"\eta = 0.2",linewidth = 2, thickness_scaling = 1)
plot!(tau_entropy_03[1:200],label=L"\eta = 0.3",linewidth = 2, thickness_scaling = 1)
plot!(tau_entropy_04[1:200],label=L"\eta = 0.4",linewidth = 2, thickness_scaling = 1)
plot!(tau_entropy_05[1:200],label=L"\eta = 0.5",linewidth = 2, thickness_scaling = 1)
plot!(tau_entropy_06[1:200],label=L"\eta = 0.6",linewidth = 2, thickness_scaling = 1)

savefig("fig4.png")

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

times = LinRange(0.05,10,TT)

for i in 1:length(times)
    t = times[i]
    for j in 1:M
        idx = findall(x->x <= t,delta_cum[:,j,1])
        if isempty(idx)
            time_values[i,j,1] = tau_entropy_values_00[1,j]
        else
            time_values[i,j,1] = tau_entropy_values_00[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,2])
        if isempty(idx)
            time_values[i,j,2] = tau_entropy_values_01[1,j]
        else
            time_values[i,j,2] = tau_entropy_values_01[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,3])
        if isempty(idx)
            time_values[i,j,3] = tau_entropy_values_02[1,j]
        else
            time_values[i,j,3] = tau_entropy_values_02[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,4])
        if isempty(idx)
            time_values[i,j,4] = tau_entropy_values_03[1,j]
        else
            time_values[i,j,4] = tau_entropy_values_03[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,5])
        if isempty(idx)
            time_values[i,j,5] = tau_entropy_values_04[1,j]
        else
            time_values[i,j,5] = tau_entropy_values_04[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,6])
        if isempty(idx)
            time_values[i,j,6] = tau_entropy_values_05[1,j]
        else
            time_values[i,j,6] = tau_entropy_values_05[idx[end],j]
        end
        idx = findall(x->x <= t,delta_cum[:,j,7])
        if isempty(idx)
            time_values[i,j,7] = tau_entropy_values_06[1,j]
        else
            time_values[i,j,7] = tau_entropy_values_06[idx[end],j]
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

default(palette = palette(:viridis,9, rev = false))
p2 = plot(times,time_mean[:,1,1],xlabel="Time [s]",ylabel="Posterior Entropy [bit]", legend = false,linewidth = 2, thickness_scaling = 1, size = (400, 300), dpi=300)
plot!(times,time_mean[:,1,2],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,3],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,4],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,5],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,6],linewidth = 2, legend = false)
plot!(times,time_mean[:,1,7],linewidth = 2, legend = false)

plot!([0,0.1,0.2,0.3,0.4,0.5,0.6],info, inset = (1, bbox(0.6,0.0,0.4,0.4)), subplot = 2, legend = false,linewidth = 2,xlabel=L"\eta",ylabel="[bit/s]",color="black")

savefig("fig5.png")

