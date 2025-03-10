@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie, LaTeXStrings"

using Printf
using CairoMakie
using LaTeXStrings
using Statistics
using Oceananigans
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days

filename = "./OUTPUTS/RB_gpu_simulation"

@info"Setting up model"

const Nx = 64     # number of points in each of horizontal directions
const Nz = 32          # number of points in the vertical direction

const Lx = 32     # (m) domain horizontal extents
const Lz = 8          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (0,Lz),
                       topology = (Bounded, Flat, Bounded)
)

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
avg_T_timeseries = FieldTimeSeries(filename * ".jld2", "avg_T")
times = T_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (3200, 1600))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)",
               aspect = DataAspect()
)

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+v^2}$", axis_kwargs...)
ax_avg_T = Axis(fig[2,3]; title = L"Average Temperature over $x$", xlabel = "T", ylabel = "z(m)")

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]
avg_T = @lift vec(dropdims(avg_T_timeseries[:, :, :, $n], dims=(1,2)))

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))

xlims!(ax_avg_T, Tlims)
ylims!(ax_avg_T, 0, Lz)
z_vec = LinRange(0, Lz, Nz)
lines!(ax_avg_T, avg_T, color=:red)

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)

w_timeseries = FieldTimeSeries(filename * ".jld2", "w")

zrange = axes(w_timeseries, 3)  # the actual valid indices in the 3rd dimension
w_center_timeseries = 0.5 .* (
    w_timeseries[:, :, zrange[1:end-1], :] .+ w_timeseries[:, :, zrange[2:end], :]
)

wT_timeseries = w_center_timeseries .* T_timeseries

wT_avg_timeseries = mean(wT_timeseries, dims=(1))

wT_avg = @lift vec(dropdims(wT_avg_timeseries[:, :, :, $n], dims=(1,2)))

ax_wT = Axis(fig[3,3];
    title = L"Vertical Heat Flux, $wT$ (averaged over $x$)",
    xlabel = "wT", ylabel = "z (m)"
)

wTlims = (minimum(wT_avg_timeseries), maximum(wT_avg_timeseries))
xlims!(ax_wT, wTlims)
ylims!(ax_wT, 0, Lz)

lines!(ax_wT, z_vec, wT_avg, color=:red)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=16) do i
    n[] = i
end
#=
#Nusselt num

avg_wT = mean(wT_timeseries)

avged_wT = Lz * avg_wT / (κ * Δ)

@info "⟨wT⟩ = $avged_wT"
=#
# Compute mean wT over x, y, z
wT_avg_timeseries_2 = mean(wT_timeseries, dims=(1,2,3))  

# Create figure for Nusselt number vs time
fig_Nu = Figure(size = (1200, 600))
ax_Nu = Axis(fig_Nu[1, 1];
    title = "Nusselt Number vs Time",
    xlabel = "Time (days)", ylabel = "Nusselt Number"
)

Nu = @lift 1 .+ vec(dropdims(wT_avg_timeseries_2[:, :, :, $n], dims=(3)))

Nulims = (minimum(wT_avg_timeseries_2)+1, maximum(wT_avg_timeseries_2)+1)
ylims!(ax_Nu, Nulims)

# Plot Nusselt number vs time
times_days = times ./ days
lines!(ax_Nu, times_days, Nu)

# Save figure
save("./OUTPUTS/Nusselt_number_vs_time.png", fig_Nu)