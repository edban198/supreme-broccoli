using Printf
using CairoMakie
using LaTeXStrings
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Units: kilometers, kilometer, meter, meters

filename = "OUTPUTS/RB_gpu_simulation"

const Lx = 200kilometers     # (m) domain horizontal extents
const Lz = 2000meters          # (m) domain depth

const refinement = 1.2 # controls spacing near surface (higher means finer spaced)
const stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(GPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (-Lz,0),
                       topology = (Periodic, Flat, Bounded)
)


@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
ω_timeseries = FieldTimeSeries(filename * ".jld2", "ω")
avg_T_timeseries = FieldTimeSeries(filename * ".jld2", "avg_T")
times = T_timeseries.times

xT, zT = nodes(T_timeseries)
xs, zx = nodes(s_timeseries)
xω, zω = nodes(ω_timeseries)

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,1200))

axis_kwargs = (xlabel = "x (km)", ylabel = "z (m)"
)

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+v^2}$", axis_kwargs...)
ax_ω = Axis(fig[4,1]; title = L"Vorticity, $\omega = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}$", axis_kwargs...)
ax_avg_T = Axis(fig[5,1]; title = L"Average Temperature over $x$", xlabel = "T", ylabel = "z(m)")

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]
ω = @lift ω_timeseries[$n]
avg_T = @lift avg_T_timeseries[$n]

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))
ωlims = (minimum(abs, interior(ω_timeseries)), maximum(abs, interior(ω_timeseries)))

# Set axis limits explicitly to match your domain
xlims!(ax_T, 0, Lx)
ylims!(ax_T, -Lz, 0)

xlims!(ax_s, 0, Lx)
ylims!(ax_s, -Lz, 0)

xlims!(ax_ω, 0, Lx)
ylims!(ax_ω, -Lz, 0)

xlims!(ax_avg_T, Tlims)

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)
hm_ω = heatmap!(ax_ω, ω; colormap = :balance, colorrange = ωlims)
Colorbar(fig[4,2])
lines!(ax_avg_T, avg_T)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=32) do i
    n[] = i
end
