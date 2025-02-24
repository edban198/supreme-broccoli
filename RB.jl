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

const Nx = 256     # number of points in each of horizontal directions
const Nz = 256          # number of points in the vertical direction

const Lx = 64     # (m) domain horizontal extents
const Lz = 64          # (m) domain depth

grid = RectilinearGrid(GPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (-Lz,0),
                       topology = (Bounded, Flat, Bounded)
)

# Buoyancy that depends on temperature:
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4), constant_salinity=0)

const Δ = 1e-3
const Γ = 1e-6
#RB1
T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = ValueBoundaryCondition(20+Δ))
#RB2
#T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Γ), bottom = FluxBoundaryCondition(Γ))
#RB3
#T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = FluxBoundaryCondition(Γ))

const g = buoyancy.gravitational_acceleration
const α = buoyancy.equation_of_state.thermal_expansion

const ν = 1e-3
const κ = 1e-6

const Ra = g * α * Δ * Lz^3 / (ν * κ)
const Pr = ν/κ

#=
const Ra = 1e12
const Pr = 1

const ν = sqrt(g * α * Δ * Lz^3 / (Pr * Ra))
const κ = sqrt(g * α * Δ * Lz^3 * Pr / Ra)
=#

@info "Ra = $Ra, Pr = $Pr, ν = $ν, κ = $κ"

closure = ScalarDiffusivity(ν=ν,κ=κ)

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
Tᵢ(x, z) = 20

# Velocity initial condition:
uᵢ(x, z) = 1e-6 * Ξ(z)
#uᵢ(x, z) = 0

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ)

# Setting up sim

simulation = Simulation(model, Δt=2seconds, stop_time = 30day)

wizard = TimeStepWizard(cfl=1.1, max_Δt=5seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

# Output

u,v,w = model.velocities

outputs = (T = model.tracers.T,
    s = sqrt(model.velocities.u^2 + model.velocities.w^2)
)

const data_interval = 2minutes

simulation.output_writers[:simple_outputs] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(data_interval),
                     filename = filename * ".jld2",
                     overwrite_existing = true
)

@info"Running the simulation..."
run!(simulation)

@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
times = T_timeseries.times

xT, zT = nodes(T_timeseries)
xs, zs = nodes(s_timeseries)

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1200,800))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", aspect = DataAspect())

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
ax_s = Axis(fig[2,3]; title = L"Speed, $s = \sqrt{u^2+w^2}$", axis_kwargs...)

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)

hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[2,4], hm_s)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, :], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=64) do i
    n[] = i
end

#Save end image

@info "Saving final frame..."
save(filename * "_final.png", fig)