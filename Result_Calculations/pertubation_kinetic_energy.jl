#here I will write the code to calculate the pertubation kinetic energy
using Oceananigans
using Statistics
using CairoMakie

filename = "/Users/loganknudsen/Documents/UMD_Research/BottomBoundaryLayer/PSI.jld2"

u_timeseries = interior(FieldTimeSeries(filename, "u"))
v_timeseries = interior(FieldTimeSeries(filename, "v"))
w_timeseries =interior(FieldTimeSeries(filename, "w"))
t_s = FieldTimeSeries(filename, "u").times

timesteps = size(u_timeseries)[4] # pulls the length of the times series wrt time

global dot_prod = [mean(u_timeseries[:,:,:,1])^2+mean(v_timeseries[:,:,:,1])^2+mean(w_timeseries[:,:,:,1])^2]

for i = 2:(timesteps)
    ui_mean = mean(u_timeseries[:,:,:,i])^2
    vi_mean = mean(v_timeseries[:,:,:,i])^2
    wi_mean = mean(w_timeseries[:,:,:,i])^2
    global dot_prod = [dot_prod; ui_mean+vi_mean+wi_mean]
end
KE = 0.5*dot_prod

fig = Figure()
axis = (xlabel="Time", ylabel = "Perturbation Kinetic Energy")

display(plot(t_s,KE; axis))
