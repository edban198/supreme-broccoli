# dipole_shallow_water_turbulence.jl

using Oceananigans
using Oceananigans.OutputWriters  # for JLD2Writer
using JLD2
using Printf
using CairoMakie

# 1. Set up the grid and model
grid = RectilinearGrid(size  = (64, 32),
                       x     = (0, 2π),
                       y     = (-10, 10),
                       topology = (Periodic, Bounded, Flat))

model = ShallowWaterModel(grid=grid,
                          gravitational_acceleration=1.0,
                          closure=nothing)

# 2. Initial dipole in u,h fields
A₀ = 1.0; α₀ = 1.0; x₀ = π; y₀ = 0.5; H = 15.0

uᵢ(x,y) = A₀*2*(y-y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
          A₀*2*(y+y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y+y₀)^2))

vᵢ(x,y) = -A₀*2*(x-x₀)*α₀*(exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
                        exp(-α₀*((x-x₀)^2 + (y+y₀)^2)))

h̄(x,y)  = H
uhᵢ(x,y) = uᵢ(x,y)*H
vhᵢ(x,y) = vᵢ(x,y)*H

set!(model, uh=uhᵢ, vh=vhᵢ, h=h̄)

# 3. Define the vorticity Field
#    ω = ∂ₓ(vh/h) - ∂ᵧ(uh/h)
vorticity = Field(∂x(model.solution.vh ./ model.solution.h) .-
                  ∂y(model.solution.uh ./ model.solution.h))

# 4. Create the simulation and attach a JLD2Writer
simulation = Simulation(model, Δt=1e-2, stop_time=10)

filename = "shallow_water_turbulence"

simulation.output_writers[:fields] = JLD2Writer(model,
                                                (; vorticity),
                                                schedule = TimeInterval(0.1),
                                                filename = filename * ".jld2",
                                                overwrite_existing = true)
# ↳ writes a dataset "vorticity" every 0.1 time units :contentReference[oaicite:2]{index=2}

@info "Starting dipole simulation..."
run!(simulation)

# 5. Load with FieldTimeSeries and animate in CairoMakie
ω_timeseries = FieldTimeSeries(filename * ".jld2", "vorticity")
times          = ω_timeseries.times

set_theme!(Theme(fontsize = 18))
fig = Figure(resolution = (1200, 800))

axis_kwargs = (xlabel = "x",
               ylabel = "y",
               limits = ((0, 2π), (-10, 10)),
               aspect = AxisAspect(1))

ax  = Axis(fig[1, 1]; title = "Vorticity", axis_kwargs...)
n   = Observable(1)
ω   = @lift ω_timeseries[$n]           # lift the nth snapshot

hm  = heatmap!(ax, ω; colormap=:balance, colorrange = (-1, 1))
Colorbar(fig[1, 2], hm)

# record an MP4
frames = 1:length(times)
@info "Recording animation..."
record(fig, filename * ".mp4", frames; framerate = 15) do i
    n[]      = i
    ax.title = @sprintf("Time = %1.2f", times[i])
end