@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie"

using Printf
using CairoMakie
using LaTeXStrings
using Statistics
using Oceananigans
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Units: kilometers, kilometer, meter, meters

filename = "OUTPUTS/cpu_wind_simulation"

@info"Setting up model"

const Nx = 32     # number of points in each of horizontal directions
const Nz = 32          # number of points in the vertical direction

const Lx = 1000meters     # (m) domain horizontal extents
const Lz = 1000meters          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (-Lz,0),
                       topology = (Periodic, Flat, Bounded)
)

# SeawaterBuoyancy:
#buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())
buoyancy = BuoyancyTracer()

closure = ScalarDiffusivity(ν=1e-6, κ=1.4e-7)

const ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
const u₁₀ = 5    # m s⁻¹, average wind velocity 10 meters above the ocean (orig 10)
const cᴰ = 2.58e-3 # dimensionless drag coefficient
const ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

const τx = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

const f = 1e-4 # s⁻¹, Coriolis parameter
const ωₜ = 0.95 * f

const k = 2π / Lx # m⁻¹, horizontal wavenumber

const tₑ = 20days
inertial_wave(x,t) = t ≤ tₑ ? τx*sin(ωₜ*t) : 0.0

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(inertial_wave))
#=
heaviside(x) = ifelse(x<0, zero(x), one(x))

const H = grid.Lz

function bottom_mask_func(z)
    sponge_one = -H/4
    sponge_zero = sponge_one + H/10
    return heaviside(-(z-sponge_zero)) * (z-sponge_zero)^2 / (sponge_one-sponge_zero)^2
end

function mask_tanh(z)
    return - (tanh((z+H/2)) - 1)
end

sponge = Relaxation(rate = 1/30minutes, mask = bottom_mask_func, target=0)
=#
model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:b),
                            boundary_conditions = (u=u_bcs,),
                            closure = closure
)

#Set buoyancy with N²
const N² = 1e-5
b₀(x, z) = N² * (z)
ξ(x, z) = exp(-(x^2  + (z + 50)^2)/20)
buoy(x,z) = b₀(x,z)
set!(model, b = buoy)

# Initial conditions - nothing at the moment

# Random noise
Ξ(z) = randn()

# Temperature initial condition:
Tᵢ(x, z) = 20
#a stable density gradient with random noise superposed - 
#Tᵢ(x, z) = 20 - dTdz * z + dTdz * model.grid.Lz * 1e-5 * Ξ(z)

Sᵢ(x, z) = 35

# Velocity initial condition:
uᵢ(x, z) = 1e-6 * Ξ(z)
#uᵢ(x, z) = 0

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ)

# Setting up sim
simulation = Simulation(model, Δt=30seconds, stop_time = 2days)

wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=30seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

# Output
u,v,w = model.velocities

outputs = (s = sqrt(model.velocities.u^2 + model.velocities.w^2),
           b′ = model.tracers.b - b₀,
)

data_interval = 10minutes

simulation.output_writers[:simple_outputs] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(data_interval),
                     filename = filename * ".jld2",
                     overwrite_existing = true
)

@info"Running simulation"
run!(simulation)

s_timeseries = FieldTimeSeries(filename * ".jld2", "s")

@info s_timeseries
#record movie
frames = 1:length(times)
@info "Making an animation..."
CairoMakie.record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end