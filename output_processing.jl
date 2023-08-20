using Oceananigans
using GLMakie#master

filename = "/Users/loganknudsen/Documents/UMD_Research/BottomBoundaryLayer/internal_wave.jld2"

set_theme!(Theme(fontsize = 24))

fig = Figure(resolution = (600, 600))

ax = Axis(fig[2, 1]; xlabel = "y", ylabel = "z",
          limits = ((0, 3000), (-200, 0)), aspect = AxisAspect(1))


# Next, we load `w` data with `FieldTimeSeries` of `w` and make contour
# plots of vertical velocity. We use Makie's `Observable` to animate the data.
# To dive into how `Observable`s work, refer to
# [Makie.jl's Documentation](https://makie.juliaplots.org/stable/documentation/nodes/index.html).

n = Observable(1)

w_timeseries = FieldTimeSeries(filename, "u")
x, y, z = nodes(w_timeseries)

w = @lift interior(w_timeseries[$n], 1, :, :)
w_lim = 1

@info "Data has been stored into w"

plt = GLMakie.heatmap!(ax, y, z, w,
          colormap = :balance,
          colorrange = (-w_lim, w_lim)
          )
Colorbar(fig[2,2],plt)
@info "Data has been mapped"

title = @lift L"$\frac{tf}{2\pi}$ = " * string(round(w_timeseries.times[$n], digits=2))
fig[1, 1] = Label(fig, title, fontsize=24, tellwidth=false)

# And, finally, we record a movie.

@info "Compiling movie"

tlength=length(w_timeseries.times)

frames = 1:tlength
record(fig, "/Users/loganknudsen/Documents/UMD_Research/BottomBoundaryLayer/internal_waves.mp4", frames, framerate=8) do j
    msg = string("Plotting frame ", j, "of ", frames[end])
    # @info msg * "\r"
    n[] = j
end


# ![](internal_wave.mp4)
