@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie"

using Printf
using CUDA
using CairoMakie
using LaTeXStrings
using Statistics
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans
using Oceananigans.Units: seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Units: kilometers, kilometer, meter, meters

filename = "OUTPUTS/cpu_wave_propagation"

@info"Setting up model"

const Nx = 512     # number of points in each of horizontal directions
const Nz = 196          # number of points in the vertical direction

const Lx = 5kilometers     # (m) domain horizontal extents
const Lz = 1000meters          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (-Lz,0),
                       topology = (Periodic, Flat, Bounded)
)

closure = ScalarDiffusivity()
#buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState())

# Revised parameters for wave generation at the surface
const C = 1e-4    # Wave amplitude [m/s]
const a = 0.01    # Inverse width of surface forcing layer [1/m]
const z₀ = -100    # Depth of forcing center [m] (50m below surface)
const λ = 100
const k_z = 2π/λ # Vertical wavenumber (50m wavelength)
const ω = 1e-3    # Angular frequency [rad/s]

# Surface-localized forcing function (1 at surface, 0 at depth)
step_func(z) = 0.5 * (1 + tanh(a * (z - z₀)))

# Divergence-free forcing for incompressible flow
u_forcing(x, z, t) =  C * cos(k_z * z + ω * t) * step_func(z)
w_forcing(x, z, t) = -C * sin(k_z * z + ω * t) * step_func(z)

# Free-slip at bottom, wave-permissive at top
u_bcs = FieldBoundaryConditions(
    top = FluxBoundaryCondition(0.0),    # Allow vertical velocity
    bottom = ValueBoundaryCondition(0.0) # No-slip at bottom
)

model = NonhydrostaticModel(; grid,
    closure = closure,
    boundary_conditions = (u=u_bcs,),
    forcing = (u=u_forcing, w=w_forcing)
)

# Velocity initial condition:
# Start from rest instead of random noise
set!(model, u=0, w=0)

# Match time-stepping to wave frequency
simulation = Simulation(model, Δt=2π/(10ω), stop_time=20days)

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
           w = model.velocities.w
)

const data_interval = 6minutes

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
times = w_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,1200))

axis_kwargs = (xlabel = "x (km)", ylabel = "z (m)"
)

ax_w = Axis(fig[2,1]; title = L"Veritcal Velocity, $w$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+w^2}$", axis_kwargs...)
ax_ω = Axis(fig[4,1]; title = L"Vorticity, $\omega = \frac{\partial u}{\partial z} - \frac{\partial w}{\partial x}$", axis_kwargs...)

n = Observable(1)

w = @lift w_timeseries[$n]
s = @lift s_timeseries[$n]
ω = @lift ω_timeseries[$n]

wlims = (minimum(interior(w_timeseries)), maximum(interior(w_timeseries)))
slims = (minimum(interior(s_timeseries)), maximum(interior(s_timeseries)))
ωlims = (minimum(interior(ω_timeseries)), maximum(interior(ω_timeseries)))

@info wlims
@info slims
@info ωlims

# Set axis limits explicitly to match your domain
xlims!(ax_w, 0, Lx)
ylims!(ax_w, -Lz, 0)

xlims!(ax_s, 0, Lx)
ylims!(ax_s, -Lz, 0)

xlims!(ax_ω, 0, Lx)
ylims!(ax_ω, -Lz, 0)

hm_w = heatmap!(ax_w, w; colormap = :balance, colorrange = wlims)
Colorbar(fig[2,2], hm_w)

hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)

hm_ω = heatmap!(ax_ω, ω; colormap = :balance, colorrange = ωlims)
Colorbar(fig[4,2], hm_ω)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
CairoMakie.record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end