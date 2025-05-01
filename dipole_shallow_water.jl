# dipole_shallow_water_with_diagnostic.jl

using Oceananigans
using Oceananigans.OutputWriters  # for JLD2Writer
using Oceananigans.Diagnostics     # for Computation & Diagnostics
using JLD2
using CairoMakie
using Printf

# 1. Grid & model
grid = RectilinearGrid(size     = (64, 32),
                       x        = (0, 2π),
                       y        = (-10, 10),
                       topology = (Periodic, Bounded, Flat))

model = ShallowWaterModel(grid=grid,
                          gravitational_acceleration=1.0,
                          closure=nothing)

# 2. Initial dipole in uh, vh, h
A₀, α₀, x₀, y₀, H = 1.0, 1.0, π, 0.5, 15.0

u₀(x,y) = A₀*2*(y-y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
         A₀*2*(y+y₀)*α₀*exp(-α₀*((x-x₀)^2 + (y+y₀)^2))

v₀(x,y) = -A₀*2*(x-x₀)*α₀*(exp(-α₀*((x-x₀)^2 + (y-y₀)^2)) -
                        exp(-α₀*((x-x₀)^2 + (y+y₀)^2)))

set!(model,
     uh = (x,y) -> u₀(x,y)*H,
     vh = (x,y) -> v₀(x,y)*H,
     h  = (x,y) -> H)

# 3. Define vorticity operation & computation
#    Velocities u,v are interpolated to cell centers automatically.
u, v = model.velocities
vorticity_op = ∂x(v) - ∂y(u)                   # Operation: ∂ₓv − ∂ᵧu :contentReference[oaicite:3]{index=3}

ω = Field(Face, Face, Cell, model.architecture, grid)  # storage Field
vorticity_comp = Computation(vorticity_op, ω)          # hook up op → ω
push!(model.diagnostics, vorticity_comp)                # run at each RK stage :contentReference[oaicite:4]{index=4}

# 4. Simulation + output writer
simulation = Simulation(model,
                        Δt       = 1e-2,
                        stop_time= 10.0)

simulation.output_writers[:vort] = JLD2Writer(model,
                                               model.diagnostics,       # dumps every diagnostic
                                               schedule             = TimeInterval(0.1),
                                               filename             = "shallow_water_dipole.jld2",
                                               overwrite_existing   = true)  # :contentReference[oaicite:5]{index=5}

@info "Running shallow-water dipole..."
run!(simulation)

# 5. Animate with CairoMakie
ts = FieldTimeSeries("shallow_water_dipole.jld2", "vorticity")
times = ts.times

fig = Figure(resolution=(1200,800))
ax  = Axis(fig[1,1],
           xlabel="x", ylabel="y",
           limits=((0,2π),(-10,10)), aspect=AxisAspect(1))

n = Observable(1)
ωsnap = @lift ts[$n]   # the nth vorticity snapshot

hm = heatmap!(ax, ωsnap; colormap=:balance, colorrange=(-1,1))
Colorbar(fig[1,2], hm)

record(fig, "dipole_vorticity.mp4", 1:length(times); framerate=15) do i
  n[]      = i
  ax.title = @sprintf("t = %1.2f", times[i])
end
