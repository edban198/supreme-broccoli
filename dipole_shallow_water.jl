using Oceananigans
using Oceananigans.OutputWriters: JLD2OutputWriter
using JLD2
using Printf

# Model setup
grid = RectilinearGrid(size=(64, 32), x=(0, 2π), y=(-10, 10), topology=(Periodic, Bounded, Flat))

model = ShallowWaterModel(; grid,
    gravitational_acceleration=1.0,
    closure=nothing
)

# Initial conditions
A₀ = 1.0
α₀ = 1.0
x₀ = π
y₀ = 0.5
H = 15.0

uᵢ(x, y) = A₀ * 2(y - y₀) * α₀ * exp(-α₀*((x - x₀)^2 + (y - y₀)^2)) - A₀ * 2(y + y₀) * α₀ * exp(-α₀*((x - x₀)^2 + (y + y₀)^2))
vᵢ(x, y) = -A₀ * 2(x - x₀) * α₀ * (exp(-α₀*((x - x₀)^2 + (y - y₀)^2)) - exp(-α₀*((x - x₀)^2 + (y + y₀)^2)))

h̄(x, y) = H
uhᵢ(x, y) = uᵢ(x, y) * h̄(x, y)
vhᵢ(x, y) = vᵢ(x, y) * h̄(x, y)

set!(model, uh=uhᵢ, vh=vhᵢ, h=h̄)

# Simulation setup
simulation = Simulation(model, Δt=1e-2, stop_time=10)

# Output setup
ω = Field(∂x(model.solution.vh/model.solution.h) - ∂y(model.solution.uh/model.solution.h))
outputs = Dict(:vorticity => ω)

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs;
    filename="dipole_shallow_water.jld2",
    schedule=TimeInterval(0.1),
    overwrite_existing=true
)

# Run simulation
@info "Starting simulation..."
run!(simulation)

using JLD2, CairoMakie

function visualize()
    file = jldopen("dipole_shallow_water.jld2")
    times = parse.(Float64, keys(file["timeseries/t"]))
    
    fig = Figure(resolution=(1200, 800))
    ax = Axis(fig[1, 1], xlabel="x", ylabel="y")
    
    n = Observable(1)
    ω = @lift file["timeseries/vorticity"][:, :, 1, $n]
    hm = heatmap!(ax, ω, colorrange=(-1, 1), colormap=:balance)
    Colorbar(fig[1, 2], hm)
    
    record(fig, "vorticity_animation.mp4", 1:length(times), framerate=24) do i
        n[] = i
    end
    
    close(file)
end

visualize()