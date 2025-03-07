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

const Nx = 128     # number of points in each of horizontal directions
const Nz = 128          # number of points in the vertical direction

const Lx = 8     # (m) domain horizontal extents
const Lz = 8          # (m) domain depth

grid = RectilinearGrid(GPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (0,Lz),
                       topology = (Bounded, Flat, Bounded)
)

# Buoyancy that depends on temperature:
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=0)

#Set values
const R = 1707.76 * 1.01
const Pr = 7.0
const ν = 1.04e-5
const κ = ν / Pr
const g = buoyancy.gravitational_acceleration
const α = buoyancy.equation_of_state.thermal_expansion
const Δ = ν * κ * R / (g * α * Lz^3)

T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0), bottom = ValueBoundaryCondition(Δ))

closure = ScalarDiffusivity(ν=ν,κ=κ)

const f = 10 * κ / Lz^2

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:T),
                            closure = closure,
                            boundary_conditions = (; T=T_bcs)
)

# Initial conditions

# Random noise
Ξ(z) = randn()

# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, z) = Δ * (1 - z/Lz)

# Velocity initial condition:
uᵢ(x, z) = 1e-6 * Ξ(z)
#uᵢ(x, z) = 0

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ)

# Setting up sim

simulation = Simulation(model, Δt=10seconds, stop_time = 30days)

wizard = TimeStepWizard(cfl=1.1, max_Δt=20seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

# Output
# Define outputs to save w and T directly
outputs = (
    w = model.velocities.w,
    T = model.tracers.T,
    avg_T = mean(model.tracers.T, dims=(1,2)),
    s = sqrt(model.velocities.u^2 + model.velocities.w^2)  # Optional
)

const data_interval = 2minutes

simulation.output_writers[:full_outputs] = JLD2OutputWriter(
    model, outputs,
    schedule = TimeInterval(data_interval),
    filename = filename * ".jld2",
    overwrite_existing = true
)

@info"Running the simulation..."
run!(simulation)
@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
avg_T_timeseries = FieldTimeSeries(filename * ".jld2", "avg_T")
times = T_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (800,1800))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)",
               aspect = DataAspect()
)

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+v^2}$", axis_kwargs...)
ax_avg_T = Axis(fig[4,1]; title = L"Average Temperature over $x$", xlabel = "T", ylabel = "z(m)")

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]
avg_T = @lift avg_T_timeseries[$n]

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))

xlims!(ax_avg_T, Tlims)

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)
lines!(ax_avg_T, avg_T)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=16) do i
    n[] = i
end
#=
@info "Plotting vertical flux"

# Get dimensions from ACTUAL data
w_sample = w_timeseries[1][:, 1, :]
actual_Nz = size(T_timeseries[1], 3)  # z-cells in T
actual_Nfaces = size(w_sample, 2)     # z-faces in w

# Initialize based on ACTUAL dimensions
wT_horiz_avg_time_accumulator = zeros(actual_Nfaces - 2)  # Interior faces

# Vertical coordinates using ACTUAL sizes
zF = znodes(grid, Face()) |> collect
z_plot = zF[2:end-1]  # Interior faces for plotting

for (n, t) in enumerate(times)
    w = w_timeseries[n][:, 1, :]        # [x, z] faces
    T_current = T_timeseries[n][:, 1, :] # [x, z] centers
    
    # Interior faces only (exact slicing adapts to actual grid)
    T_at_w_faces = 0.5*(T_current[:, 1:end-1] + T_current[:, 2:end])
    w_interior = w[:, 2:end-1]
    
    # Dimension safeguard
    @assert size(T_at_w_faces) == size(w_interior) "Mismatched flux dimensions"
    
    wT = w_interior .* T_at_w_faces
    wT_horiz_avg = mean(wT, dims=1)
    
    wT_horiz_avg_time_accumulator .+= vec(wT_horiz_avg)
end

wT_time_avg = wT_horiz_avg_time_accumulator / length(times)

# Plot with dynamically determined z-coordinates
fig = Figure()
ax = Axis(fig[1,1], xlabel="⟨wT⟩", ylabel="z (m)")
lines!(ax, wT_time_avg, z_plot)
savefig(filename * "_flux.png")
=#