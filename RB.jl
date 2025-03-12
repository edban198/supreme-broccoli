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

const Nx = 32     # number of points in each of horizontal directions
const Nz = 16          # number of points in the vertical direction

const Lx = 32     # (m) domain horizontal extents
const Lz = 8          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (0,Lz),
                       topology = (Periodic, Flat, Bounded)
)

# Buoyancy that depends on temperature:
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=0)

#Set values
const R = 657.5 * 2
const Pr = 7.0
const ν = 1.04e-4
const κ = ν / Pr
const g = buoyancy.gravitational_acceleration
const α = buoyancy.equation_of_state.thermal_expansion
const Δ = ν * κ * R / (g * α * Lz^3) 
t_ff = sqrt(Lz / (g * α * Δ))
t_ff_days = t_ff / (3600 * 24)
@info "Freefall time in days ~ $t_ff_days"

T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0), bottom = ValueBoundaryCondition(Δ))

τx = -1e-6
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τx))

closure = ScalarDiffusivity(ν=ν,κ=κ)

const f = 10 * κ / Lz^2

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=5),
                            tracers = (:T),
                            closure = closure,
                            boundary_conditions = (; T=T_bcs, u_bcs)
)

# Initial conditions

# Random noise
Ξ(x,z) = randn()

# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, z) = Δ * (1 - z/Lz)

# Velocity initial condition:
uᵢ(x, z) = 1e-6 * Ξ(x,z)
#uᵢ(x, z) = 0

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ)

# Setting up sim

simulation = Simulation(model, Δt=10minutes, stop_time = 10days)

wizard = TimeStepWizard(cfl=1.1, max_Δt=1hour)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(100))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))
#=
@info "running sim..."
run!(simulation)
=#

# OUTPUTS
outputs = (
    w = model.velocities.w,
    T = model.tracers.T,
    #avg_T = mean(model.tracers.T, dims=(1,2)),
    s = sqrt(model.velocities.u^2 + model.velocities.w^2)
)

const data_interval = 2minutes

simulation.output_writers[:full_outputs] = JLD2OutputWriter(
    model, outputs,
    schedule = TimeInterval(data_interval),
    filename = filename * ".jld2",
    overwrite_existing = true
)
#=
@info"Restarting the simulation..."
simulation.stop_time = 80days
=#
run!(simulation)
@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
#avg_T_timeseries = FieldTimeSeries(filename * ".jld2", "avg_T")
times = T_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1600,1000))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)",
               aspect = DataAspect()
)

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+v^2}$", axis_kwargs...)
#ax_avg_T = Axis(fig[2,3]; title = L"Average Temperature over $x$", xlabel = "T", ylabel = "z(m)")

n = Observable(1)

T = @lift T_timeseries[$n]
s = @lift s_timeseries[$n]
#avg_T = @lift vec(dropdims(avg_T_timeseries[:, :, :, $n], dims=(1,2)))

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))
#=
xlims!(ax_avg_T, Tlims)
ylims!(ax_avg_T, 0, Lz)
z_vec = LinRange(0, Lz, Nz)
lines!(ax_avg_Tavg_T, color=:red)
=#
hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
Colorbar(fig[3,2], hm_s)
#=
using Interpolations

zrange = axes(w_timeseries, 3)  # the actual valid indices in the 3rd dimension
w_center_timeseries = 0.5 .* (
    w_timeseries[:, :, zrange[1:end-1], :] .+ w_timeseries[:, :, zrange[2:end], :]
)

wT_timeseries = w_center_timeseries .* T_timeseries

wT_avg_timeseries = mean(wT_timeseries, dims=(1))

wT_avg = @lift vec(dropdims(wT_avg_timeseries[:, :, :, $n], dims=(1,2)))

ax_wT = Axis(fig[3,3];
    title = L"Vertical Heat Flux, $wT$ (averaged over $x$)",
    xlabel = "wT", ylabel = "z (m)"
)

wTlims = (minimum(wT_avg_timeseries), maximum(wT_avg_timeseries))
xlims!(ax_wT, wTlims)
ylims!(ax_wT, 0, Lz)

lines!(ax_wT, wT_avg, color=:red)
=#
@info "calculating Nusselt number"
#Nusselt num
w_timeseries = FieldTimeSeries(filename * ".jld2", "w")
zrange = axes(w_timeseries, 3)  # the actual valid indices in the 3rd dimension
w_center_timeseries = 0.5 .* (
    w_timeseries[:, :, zrange[1:end-1], :] .+ w_timeseries[:, :, zrange[2:end], :]
)

wT_timeseries = w_center_timeseries .* T_timeseries

avg_wT = mean(wT_timeseries)

Nu = 1 + avg_wT

@info "Nu = $Nu"

title = @lift "t = " * prettytime(times[$n]) * ", Nu = " * string(round(Nu, digits=3))
Label(fig[1, :], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=16) do i
    n[] = i
end

#=
# Compute mean wT over x, y, z
wT_avg_timeseries_2 = mean(wT_timeseries, dims=(1,2,3))  

# Create figure for Nusselt number vs time
fig_Nu = Figure(size = (1200, 600))
ax_Nu = Axis(fig_Nu[1, 1];
    title = "Nusselt Number vs Time",
    xlabel = "Time (days)", ylabel = "Nusselt Number"
)

Nu = @lift 1 .+ Lz .* vec(dropdims(wT_avg_timeseries_2[:, :, :, $n], dims=(3))) ./ (κ * Δ)

Nulims = (minimum(wT_avg_timeseries_2)+1, maximum(wT_avg_timeseries_2)+1)
ylims!(ax_Nu, Nulims)

# Plot Nusselt number vs time
times_days = times ./ days
lines!(ax_Nu, times_days, Nu)

# Save figure
save("./OUTPUTS/Nusselt_number_vs_time.png", fig_Nu)
=#