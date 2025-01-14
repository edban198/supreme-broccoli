@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie"

using Printf
using CairoMakie
using LaTeXStrings
using Statistics
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Units: kilometers, kilometer, meter, meters

filename = "OUTPUTS/RB_gpu_simulation"

@info"Setting up model"

const Nx = 512     # number of points in each of horizontal directions
const Nz = 196          # number of points in the vertical direction

const Lx = 10kilometers     # (m) domain horizontal extents
const Lz = 1000meters          # (m) domain depth

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

# Buoyancy that depends on temperature:

buoyancy = SeawaterBuoyancy(equation_of_state = TEOS10EquationOfState())

const Δ = 1e-3
const Γ = 1e-6
#RB1
#T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = ValueBoundaryCondition(20+Δ))
#RB2
#T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Γ), bottom = FluxBoundaryCondition(Γ))
#RB3
#T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(20), bottom = FluxBoundaryCondition(Γ))

#const g = buoyancy.gravitational_acceleration
#const α = buoyancy.equation_of_state.thermal_expansion

#=
const ν = 1e-3
const κ = 1e-6

const Ra = g * α * Δ * Lz^3 / (ν * κ)
const Pr = ν/κ
=#

#const Ra = 1e12
#const Pr = 1

#const ν = sqrt(g * α * Δ * Lz^3 / (Pr * Ra))
#const κ = sqrt(g * α * Δ * Lz^3 * Pr / Ra)

closure = ScalarDiffusivity()

sim_length = 50days
Δt = 20seconds

const τx = 1e-6 #wind flux

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τx))

heaviside(x) = ifelse(x<0, zero(x), one(x))

const H = grid.Lz

sponge_one = minimum(Oceananigans.Grids.znodes(grid, Face()))
sponge_zero = sponge_one + grid.Lz/10

function bottom_mask_func(x,z)
    sponge_one = -H/2
    sponge_zero = sponge_one + H/10
    return heaviside(-(z-sponge_zero)) * (z-sponge_zero)^2 / (sponge_one-sponge_zero)^2
end

sponge = Relaxation(rate = 1/30minutes, mask = bottom_mask_func, target=0)

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:b,:T,:S),
                            closure = closure,
                            boundary_conditions = (u=u_bcs,),
                            forcing = (w=sponge,)
)

# Initial conditions

# Random noise
Ξ(z) = randn()

# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, z) = 20
#Tᵢ(x, z) = 20 - dTdz * z + dTdz * model.grid.Lz * 1e-5 * Ξ(z)

# Velocity initial condition:
uᵢ(x, z) = 1e-6 * Ξ(z)
#uᵢ(x, z) = 0

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ)

# Setting up sim

simulation = Simulation(model, Δt=Δt, stop_time = sim_length)

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

u,v,w = model.velocities

outputs = (s = sqrt(model.velocities.u^2 + model.velocities.w^2),
           ω = Field(∂z(model.velocities.u) - ∂x(model.velocities.w)),
           w = model.velocities.w,
           b = model.tracers.b
)

const data_interval = 10minutes

simulation.output_writers[:simple_outputs] =
    JLD2OutputWriter(model, outputs,
                     schedule = TimeInterval(data_interval),
                     filename = filename * ".jld2",
                     overwrite_existing = true
)

@info"Running the simulation..."
run!(simulation)

@info"Plotting animation"

w_timeseries = FieldTimeSeries(filename * ".jld2", "w")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
ω_timeseries = FieldTimeSeries(filename * ".jld2", "ω")
b_timeseries = FieldTimeSeries(filename * ".jld2", "b")
times = w_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,1200))

axis_kwargs = (xlabel = "x (km)", ylabel = "z (m)"
)

ax_w = Axis(fig[2,1]; title = L"Veritcal Velocity, $w$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+w^2}$", axis_kwargs...)
ax_ω = Axis(fig[4,1]; title = L"Vorticity, $\omega = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}$", axis_kwargs...)
ax_b = Axis(fig[5,1]; title = L"Buoyancy, $b$", axis_kwargs...)

n = Observable(1)

w = @lift w_timeseries[$n]
s = @lift s_timeseries[$n]
ω = @lift ω_timeseries[$n]
b = @lift b_timeseries[$n]

wlims = (minimum(abs, interior(w_timeseries)), maximum(abs, interior(w_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))
ωlims = (minimum(abs, interior(ω_timeseries)), maximum(abs, interior(ω_timeseries)))
blims = (minimum(abs, interior(b_timeseries)), maximum(abs, interior(b_timeseries)))

@info wlims
@info slims
@info ωlims
@info blims

# Set axis limits explicitly to match your domain
xlims!(ax_w, 0, Lx)
ylims!(ax_w, -Lz, 0)

xlims!(ax_s, 0, Lx)
ylims!(ax_s, -Lz, 0)

xlims!(ax_ω, 0, Lx)
ylims!(ax_ω, -Lz, 0)

xlims!(ax_b, 0, Lx)
ylims!(ax_b, -Lz, 0)

hm_w = heatmap!(ax_w, w; colormap = :balance, colorrange = wlims)
Colorbar(fig[2,2], hm_w)

hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)

hm_ω = heatmap!(ax_ω, ω; colormap = :balance, colorrange = ωlims)
Colorbar(fig[4,2])

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end