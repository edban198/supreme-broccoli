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

const Nx = 256     # number of points in each of horizontal directions
const Nz = 128          # number of points in the vertical direction

const Lx = 5kilometers     # (m) domain horizontal extents
const Lz = 1000meters          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (-Lz,0),
                       topology = (Periodic, Flat, Bounded)
)

# SeawaterBuoyancy:
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState())

closure = ScalarDiffusivity(ν=1e-6, κ=1.4e-7)

const ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
const u₁₀ = 5    # m s⁻¹, average wind velocity 10 meters above the ocean (orig 10)
const cᴰ = 2.58e-3 # dimensionless drag coefficient
const ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

const τx = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) # m² s⁻²

const f = 1e-4 # s⁻¹, Coriolis parameter

const k = 2π / Lx # m⁻¹, horizontal wavenumber

inertial_wave(x,t) = τx

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
model = NonhydrostaticModel(; grid,
                            advection = UpwindBiased(order=5),
                            tracers = (:T,:S),
                            coriolis = FPlane(f),
                            closure = closure,
                            boundary_conditions = (u=u_bcs,)
)

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
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ, S=Sᵢ)

# Setting up sim
simulation = Simulation(model, Δt=30seconds, stop_time = 20days)

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
           ω = ∂z(model.velocities.u) - ∂x(model.velocities.w),
           w = model.velocities.w,
           u = model.velocities.u
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
u_timeseries = FieldTimeSeries(filename * ".jld2", "u")
times = w_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,1500))

axis_kwargs = (xlabel = "x (km)", ylabel = "z (m)"
)

ax_w = Axis(fig[2,1]; title = L"Veritcal Velocity, $w$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+w^2}$", axis_kwargs...)
ax_ω = Axis(fig[4,1]; title = L"Vorticity, $\omega = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}$", axis_kwargs...)
ax_u = Axis(fig[5,1]; title = L"Horizontal Velocity, $u$", axis_kwargs...)

n = Observable(1)

w = @lift w_timeseries[$n]
s = @lift s_timeseries[$n]
ω = @lift ω_timeseries[$n]
u = @lift u_timeseries[$n]

wlims = (minimum(interior(w_timeseries)), maximum(interior(w_timeseries)))
slims = (minimum(interior(s_timeseries)), maximum(interior(s_timeseries)))
ωlims = (minimum(interior(ω_timeseries)), maximum(interior(ω_timeseries)))
ulims = (minimum(interior(u_timeseries)), maximum(interior(u_timeseries)))

@info wlims
@info slims
@info ωlims
@info ulims

# Set axis limits explicitly to match domain
xlims!(ax_w, 0, Lx)
ylims!(ax_w, -Lz, 0)

xlims!(ax_s, 0, Lx)
ylims!(ax_s, -Lz, 0)

xlims!(ax_ω, 0, Lx)
ylims!(ax_ω, -Lz, 0)

xlims!(ax_u, 0, Lx)
ylims!(ax_u, -Lz, 0)

hm_w = heatmap!(ax_w, w; colormap = :speed, colorrange = wlims)
Colorbar(fig[2,2], hm_w)

hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)

hm_ω = heatmap!(ax_ω, ω; colormap = :balance, colorrange = ωlims)
Colorbar(fig[4,2], hm_ω)

hm_u = heatmap!(ax_u, u; colormap = :speed, colorrange = ulims)
Colorbar(fig[5,2], hm_u)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
CairoMakie.record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end