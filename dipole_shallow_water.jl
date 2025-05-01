# dipole_shallow_water_oceananigans_style.jl

using Oceananigans                             # loads all core Oceananigans types :contentReference[oaicite:3]{index=3}
using Oceananigans.OutputWriters              # brings in JLD2Writer, TimeInterval, etc.
using Oceananigans.Models: ShallowWaterModel
using CairoMakie
using Printf

# 1. Grid & Model Setup
grid = RectilinearGrid(size     = (64, 32),
                       x        = (0, 2π),
                       y        = (-10, 10),
                       topology = (Periodic, Bounded, Flat))   # 2D domain w/ flat y‐axis :contentReference[oaicite:4]{index=4}

model = ShallowWaterModel(grid=grid,
                          gravitational_acceleration = 1.0,
                          closure                   = nothing)

# 2. Dipole Initial Condition
A₀, α₀, x₀, y₀, H = 1.0, 1.0, π, 0.5, 15.0

u₀(x,y) = A₀*2*(y-y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
         A₀*2*(y+y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y+y₀)^2))

v₀(x,y) = -A₀*2*(x-x₀)*α₀*(exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
                        exp(-α₀*((x-x₀)^2 + (y+y₀)^2)))

set!(model,
     uh = (x, y) -> u₀(x, y)*H,
     vh = (x, y) -> v₀(x, y)*H,
     h  = (x, y) -> H)

# 3. Define Vorticity Operation
u, v = model.velocities                              # velocity Fields in vector-invariant form
ζ = ∂x(v) - ∂y(u)                                     # relative vorticity operator :contentReference[oaicite:5]{index=5}

# 4. Simulation & Output Writer
simulation = Simulation(model, Δt=1e-2, stop_time=10)

simulation.output_writers[:vorticity] = JLD2Writer(model,
                                                    (; ζ),
                                                    schedule           = TimeInterval(0.1),
                                                    filename           = "dipole_vorticity.jld2",
                                                    overwrite_existing = true)  # saves ζ every 0.1 units :contentReference[oaicite:6]{index=6}

@info "Running shallow‐water dipole simulation..."
run!(simulation)

# 5. Visualization & Animation
ζts   = FieldTimeSeries("dipole_vorticity.jld2", "ζ")  # lazy time‐series reader :contentReference[oaicite:7]{index=7}
times = ζts.times

fig = Figure(resolution = (1200, 800))
ax  = Axis(fig[1, 1];
           xlabel = "x", ylabel = "y",
           limits = ((0, 2π), (-10, 10)),
           aspect = AxisAspect(1))

n     = Observable(1)
ζsnap = @lift ζts[$n]                                   # Makie Observable + @lift :contentReference[oaicite:8]{index=8}

hm = heatmap!(ax, ζsnap; colormap=:balance, colorrange=(-1, 1))
Colorbar(fig[1, 2], hm)

record(fig, "dipole_vorticity.mp4", 1:length(times); framerate=15) do i
  n[]      = i
  ax.title = @sprintf("t = %1.2f", times[i])
end