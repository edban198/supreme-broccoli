@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie"

using Printf
using CairoMakie
using LaTeXStrings
using Statistics
using Oceananigans
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Units: kilometers, meters

filename = "convection"
filename_v = filename * "_veolcity"

@info"Setting up model"

Nx = 256     # number of points in each of horizontal directions
Nz = 96          # number of points in the vertical direction
#Nx, Nz = 64,32 #for testing quickley

Lx = 128     # (m) domain horizontal extents
Lz = 32          # (m) domain depth

refinement = 1.2 # controls spacing near surface (higher means finer spaced)
stretching = 12  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(size = (Nx, Nz),
                       x = (0,Lx),
                       z = z_faces,
                       topology = (Periodic, Flat, Bounded)
)

# Buoyancy that depends on temperature:
buoyancy = SeawaterBuoyancy(constant_salinity = 35)

dTdz = 0.01 # K m⁻¹

Δ = 1e-3
Γ = 1e-6
#RB1
T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = ValueBoundaryCondition(20+Δ))
#RB2
#T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Γ), bottom = FluxBoundaryCondition(Γ))
#RB3
#T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = FluxBoundaryCondition(Γ))

ν = 1e-3
κ = 1e-6

g = buoyancy.gravitational_acceleration
α = buoyancy.equation_of_state.thermal_expansion

Ra = g * α * Δ * Lz^3 / (ν * κ)
Pr = ν/κ 

closure_1 = AnisotropicMinimumDissipation()
closure_2 = (HorizontalScalarDiffusivity(ν=ν,κ=κ),
             VerticalScalarDiffusivity(ν=ν,κ=κ)
)
closure_3 = ScalarDiffusivity(ν=ν,κ=κ)

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:T),
                            closure = closure_3,
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

simulation = Simulation(model, Δt=20seconds, stop_time = 10days)

wizard = TimeStepWizard(cfl=1.0, max_Δt=30seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

# Output

simulation.output_writers[:temperature] =
    JLD2OutputWriter(model, model.tracers,
                     filename = filename * ".jld2",
                     schedule = TimeInterval(1minute),
                     overwrite_existing = true
)

#
u,v,w = model.velocities
s = sqrt(u^2 + w^2)

simulation.output_writers[:fields] = 
    JLD2OutputWriter(model, (; s),
                     schedule = TimeInterval(1minute),
                     filename = filename_v * ".jld2",
                     overwrite_existing = true
)

@info"Running the simulation..."
run!(simulation)

@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename_v * ".jld2", "s")
times = T_timeseries.times

xT, zT = nodes(T_timeseries)
xs, zx = nodes(s_timeseries)

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1200,800))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)",
               aspect = DataAspect()
)

ax_T = Axis(fig[2,1]; title = "Temperature", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $\sqrt{u^2+v^2}$", axis_kwargs...)

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=32) do i
    n[] = i
end