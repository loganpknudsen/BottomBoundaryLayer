using Oceananigans
using GLMakie#master

coriolis = FPlane(rotation_rate=7.292115e-5, latitude=45)
ps = (Nₒ = 81.7*coriolis.f, S = 7.7*coriolis.f, γ = 0.6, ϕ = 0, f = coriolis.f)

set_theme!(Theme(fontsize = 24))

fig = Figure(resolution = (600, 600))

ax = Axis(fig[2, 1]; xlabel = "y", ylabel = "z",
          limits = ((0, 3000), (-200, 0)), aspect = AxisAspect(1))


# Next, we load `w` data with `FieldTimeSeries` of `w` and make contour
# plots of vertical velocity. We use Makie's `Observable` to animate the data.
# To dive into how `Observable`s work, refer to
# [Makie.jl's Documentation](https://makie.juliaplots.org/stable/documentation/nodes/index.html).

n = Observable(1)
# filename = "PSI_outputs.jld2"
filename = "/Users/loganknudsen/Documents/UMD_Research/PSI.jld2"
w_timeseries = FieldTimeSeries(filename, "v")

x, y, z = nodes(w_timeseries)

t=w_timeseries.times
B_func(x, y, z, t, ps) = (ps.Nₒ^2-ps.γ*(ps.S^4/ps.f^2)*(cos(ps.ϕ)-cos(ps.f*t-ps.ϕ)))*z - ps.S^2*y #multiply by z since we integrate N^2 w.r.t z
b = [B_func(xi,yi,zi,ti,ps) for xi in x, yi in y, zi in z, ti in t]
w = @lift interior(w_timeseries[$n], 1, :, :)
w_lim = 0.001

@info "Data has been stored into w"

plt = GLMakie.heatmap!(ax, y, z, w,
          colormap = :balance,
          colorrange = (-w_lim, w_lim)
          )
# GLMakie.contour!(ax,y,z,b)
Colorbar(fig[2,2],plt)
@info "Data has been mapped"

title = @lift L"$\frac{tf}{2\pi}$ = " * string(round(w_timeseries.times[$n], digits=2))
fig[1, 1] = Label(fig, title, fontsize=24, tellwidth=false)

# And, finally, we record a movie.

@info "Compiling movie"

tlength = length(t)
frames = 1:tlength
record(fig, "/Users/loganknudsen/Documents/UMD_Research/cross_front_visualization2.mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, "of ", frames[end])
    # @info msg * "\r"
    n[] = j
end


# ![](internal_wave.mp4)
